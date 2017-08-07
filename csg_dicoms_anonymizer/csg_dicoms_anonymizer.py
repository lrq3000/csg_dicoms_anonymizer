# coding: utf-8
#
# CSG Dicoms Anonymizer
# By Stephen Larroque @ Coma Science Group, GIGA Research, University of Liege
# Creation date: 2017-02-07
# License: MIT
#
# INSTALL NOTE:
# Tested on Python 2.7.11
#
# TODO:
# * unify dicom names (if not already done)
# * unify demographics names (if not already done)
# * check if recursion ok (to anonymize MRI & PET at the same time for example).
# * convert cells to functions
# * put in a python script and use gooey (except if --cmd passed as argument)
# * freeze using pyinstaller (use @pyinstaller --noconsole csg_dicoms_anonymizer.py)
# * make a nice progress bar in gooey? add support in tqdm?
#


# In[ ]:

# Forcefully autoreload all python modules
#get_ipython().magic(u'load_ext autoreload')
#get_ipython().magic(u'autoreload 2')

from __future__ import print_function

from _version import __version__

__all__ = ['main']


# # AUX FUNCS

# In[ ]:

# Auxiliary libraries and necessary functions

import os
import re
import shutil
import sys

import zipfile
import cStringIO
import csv
import hashlib

cur_path = os.path.realpath('.')
sys.path.append(os.path.join(cur_path, 'csg_fileutil_libs'))  # for pydicom, because it does not support relative paths (yet?)

from tempfile import mkdtemp, mkstemp

import csg_fileutil_libs.pydicom as dicom
from csg_fileutil_libs.pydicom.filereader import InvalidDicomError
from csg_fileutil_libs.distance import distance
from csg_fileutil_libs.tee import Tee
from csg_fileutil_libs import argparse

from csg_fileutil_libs.aux_funcs import recwalk, replace_buggy_accents, _unidecode, _tqdm, cleanup_name, save_dict_as_csv, distance_jaccard_words_split



#***********************************
#                       AUX
#***********************************

def is_file(dirname):
    """Checks if a path is an actual file that exists"""
    if not os.path.isfile(dirname):
        msg = "{0} is not an existing file".format(dirname)
        raise ArgumentTypeError(msg)
    else:
        return dirname

def is_dir(dirname):
    """Checks if a path is an actual directory that exists"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise ArgumentTypeError(msg)
    else:
        return dirname

def is_dir_or_file(dirname):
    """Checks if a path is an actual directory that exists or a file"""
    if not os.path.isdir(dirname) and not os.path.isfile(dirname):
        msg = "{0} is not a directory nor a file".format(dirname)
        raise ArgumentTypeError(msg)
    else:
        return dirname

def get_fullpath(relpath):
    """Relative path to absolute"""
    if (type(relpath) is object or hasattr(relpath, 'read')): # relpath is either an object or file-like, try to get its name
        relpath = relpath.name
    return os.path.abspath(os.path.expanduser(relpath))

def disambiguate_names(L, dist_threshold=0.2, verbose=False):
    '''Disambiguate names in a list (ie, find all duplicate names with switched words or typos, and fix them and add them to an "alt_names" field)
    Input: list of names or list of dicts with "name" field. Output: list of dict with fields "name" and "alt_names". Alt names can then be used to do a mapping.'''
    # It's a list of dict (straight from a csv.DictReader)
    if isinstance(L[0], dict):
        res = list(L)  # copy the list of dicts (and thus all fields)
        vals = [c['name'] for c in L]  # extract names
    else:  # Convert input list to a dict
        res = [{'name': name} for name in L]
        vals = L
    for idx, c in _tqdm(enumerate(vals), total=len(vals), desc='DISAMB', unit='names', file=sys.stdout):
        for idx2, c2 in enumerate(vals[idx+1:]):
            #print(c, c2)
            if c != c2 and \
            (distance.nlevenshtein(c, c2, method=1) <= dist_threshold or distance_jaccard_words_split(c2, c, partial=True, norm=True, dist=dist_threshold) <= dist_threshold): # use shortest distance with normalized levenshtein
                if verbose:
                    print(c, c2, distance.nlevenshtein(c, c2, method=1))
                # Replace the name of the second entry with the name of the first entry
                res[idx+idx2+1]['name'] = c
                # Add the other name as an alternative name, just in case we did a mistake for example
                res[idx+idx2+1]['alt_names'] = res[idx]['alt_names'] + '/' + c2 if 'alt_names' in res[idx] else c2
    return res

def get_list_of_folders(rootpath):
    return [item for item in os.listdir(rootpath) if os.path.isdir(os.path.join(rootpath, item))]

def get_list_of_zip(rootpath):
    return [item for item in os.listdir(rootpath) if os.path.isfile(os.path.join(rootpath, item)) and item.endswith('.zip')]

def get_dcm_names_from_dir(rootpath, dcm_subj_list=None, folder_to_name=None, verbose=False):
    if dcm_subj_list is None:
        dcm_subj_list = []  # store list of subjects names from dicom files (useful for csv filtering)
    if folder_to_name is None:
        folder_to_name = {}  # store the name of the patient stored in each root folder (useful for anonymization later on)
    for subject in get_list_of_folders(rootpath):
        if verbose:
            print('- Processing subject %s' % unicode(subject, 'latin1'))
        fullpath = os.path.join(rootpath, subject)
        if not isinstance(fullpath, unicode):
            fullpath = unicode(fullpath, 'latin1')
        pts_name = None
        for dirpath, filename in recwalk(fullpath, filetype=['.dcm', '']):
            try:
                #print('* Try to read fields from dicom file: %s' % os.path.join(dirpath, filename))
                dcmdata = dicom.read_file(os.path.join(dirpath, filename), stop_before_pixels=True)  # stop_before_pixels allow for faster processing since we do not read the full dicom data, and here we can use it because we do not modify the dicom, we only read it to extract the dicom patient name
                #print(dcmdata.PatientName)
                pts_name = cleanup_name(dcmdata.PatientName)
                dcm_subj_list.append( pts_name )
                break
            except (InvalidDicomError, AttributeError) as exc:
                pass
        folder_to_name[subject] = pts_name
    return dcm_subj_list, folder_to_name

def get_dcm_names_from_zip(rootpath, dcm_subj_list=None, folder_to_name=None, verbose=False):
    if dcm_subj_list is None:
        dcm_subj_list = []  # store list of subjects names from dicom files (useful for csv filtering)
    if folder_to_name is None:
        folder_to_name = {}  # store the name of the patient stored in each root folder (useful for anonymization later on)
    # Create a temporary file (to extract one dicom from zip)
    dcmfilefh, dcmfilepath = mkstemp(suffix='.dcm')  # tesseract < 3.03 do not support "stdout" argument, so need to save into a file
    os.close(dcmfilefh)  # close file to allow writing after
    #dcmfilepath = 'tempdicomextract/tempdicom.dcm'
    #try:
    #    os.makedirs(os.path.dirname(dcmfilepath))
    #except OSError as exc:
    #    pass

    # Extract names from zipped dicom files (extract the first dicom file we can read and use its fields)
    for zipfilename in get_list_of_zip(rootpath):
        zfilepath = os.path.join(rootpath, zipfilename)
        if verbose:
            print('- Processing file %s' % zipfilename)
        with zipfile.ZipFile(zfilepath, 'r') as zipfh:
            # Extract only files, not directories (end with '/', this is standard detection in zipfile)
            zfolder = (item for item in zipfh.namelist() if item.endswith('/'))
            zfiles = (item for item in zipfh.namelist() if not item.endswith('/'))
            # Get first top folder inside zip to extract folder name (because when we will extract the zip, we need the folder name)
            try:
                folder_name = zfolder.next().strip('/')
            except StopIteration:
                folder_name = re.search('^([^\\/]+)[\\/]', zipfh.namelist()[0]).group(1)
            # Get first dicom file we can find
            pts_name = None
            for zf in zfiles:
                # Need to extract because pydicom does not support not having seek() (and zipfile in-memory does not provide seek())
                z = zipfh.read(zf) # do not use .extract(), the path can be anything and it does not support unicode (so it can easily extract to the root instead of target folder!)
                with open(dcmfilepath, 'wb') as dcmf:
                    dcmf.write(z)
                # Try to open the extracted dicom
                try:
                    if verbose:
                        print('Try to decode dicom fields with file %s' % zf)
                    dcmdata = dicom.read_file(dcmfilepath, stop_before_pixels=True)
                    pts_name = cleanup_name(dcmdata.PatientName)
                    dcm_subj_list.append( pts_name )
                    os.remove(dcmfilepath)
                    break
                except (InvalidDicomError, AttributeError) as exc:
                    continue
                except IOError as exc:
                    if 'no tag to read' in str(exc).lower():
                        continue
                    else:
                        raise
            # Add to the folder name -> dicom patient name mapping
            folder_to_name[folder_name] = pts_name
    return dcm_subj_list, folder_to_name

def dist_matrix(list1, list2, dist_threshold=0.2):
    '''Find all similar items in two lists that are below a specified distance threshold (using both letters- and words- levenshtein distances)'''
    dist_matches = {}
    for subj in list1:
        found = False
        for c in list2:
            if distance.nlevenshtein(subj, c, method=1) <= dist_threshold or distance_jaccard_words_split(subj, c, partial=True, norm=True, dist=dist_threshold) <= dist_threshold: # use shortest distance with normalized levenshtein
                if subj not in dist_matches:
                    dist_matches[subj] = []
                dist_matches[subj].append(c)
                found = True
        if not found:
            dist_matches[subj] = None
    return dist_matches


# In[ ]:

def zipwalk(zfilename):
    """Zip file tree generator.

    For each file entry in a zip archive, this yields
    a two tuple of the zip information and the data
    of the file as a StringIO object.

    zipinfo, filedata

    zipinfo is an instance of zipfile.ZipInfo class
    which gives information of the file contained
    in the zip archive. filedata is a StringIO instance
    representing the actual file data.

    If the file again a zip file, the generator extracts
    the contents of the zip file and walks them.

    Inspired by os.walk .
    Source: by Anand http://code.activestate.com/recipes/425840-zip-walker-zip-file-tree-generator/
    """

    tempdir=os.environ.get('TEMP',os.environ.get('TMP',os.environ.get('TMPDIR','/tmp')))
    
    try:
        z=zipfile.ZipFile(zfilename,'r')
        for info in z.infolist():
            fname = info.filename
            data = z.read(fname)
            extn = (os.path.splitext(fname)[1]).lower()

            if extn=='.zip':
                checkz=False
                
                tmpfpath = os.path.join(tempdir,os.path.basename(fname))
                try:
                    open(tmpfpath,'w+b').write(data)
                except (IOError, OSError),e:
                    print(e)

                if zipfile.is_zipfile(tmpfpath):
                    checkz=True

                if checkz:
                    try:
                        for x in zipwalk(tmpfpath):
                            yield x
                    except Exception, e:
                        raise
                    
                try:
                    os.remove(tmpfpath)
                except:
                    pass
            else:
                yield (info, cStringIO.StringIO(data))
    except RuntimeError, e:
        print('Runtime Error')
    except zipfile.error, e:
        raise


#***********************************
#        GUI AUX FUNCTIONS
#***********************************

# Try to import Gooey for GUI display, but manage exception so that we replace the Gooey decorator by a dummy function that will just return the main function as-is, thus keeping the compatibility with command-line usage
try:  # pragma: no cover
    import gooey
except ImportError as exc:
    # Define a dummy replacement function for Gooey to stay compatible with command-line usage
    class gooey(object):  # pragma: no cover
        def Gooey(func):
            return func
    # If --gui was specified, then there's a problem
    if len(sys.argv) > 1 and 'cmd' not in sys.argv:  # pragma: no cover
        print('ERROR: --gui specified but an error happened with lib/gooey, cannot load the GUI (however you can still use this script in commandline). Check that lib/gooey exists and that you have wxpython installed. Here is the error: ')
        raise(exc)

def conditional_decorator(flag, dec):  # pragma: no cover
    def decorate(fn):
        if flag:
            return dec(fn)
        else:
            return fn
    return decorate

def check_gui_arg():  # pragma: no cover
    """Check that the --gui argument was passed, and if true, we remove the --gui option and replace by --gui_launched so that Gooey does not loop infinitely"""
    if len(sys.argv) > 1 or '--cmd' in sys.argv:
        # DEPRECATED since Gooey automatically supply a --ignore-gooey argument when calling back the script for processing
        #sys.argv[1] = '--gui_launched' # CRITICAL: need to remove/replace the --gui argument, else it will stay in memory and when Gooey will call the script again, it will be stuck in an infinite loop calling back and forth between this script and Gooey. Thus, we need to remove this argument, but we also need to be aware that Gooey was called so that we can call gooey.GooeyParser() instead of argparse.ArgumentParser() (for better fields management like checkboxes for boolean arguments). To solve both issues, we replace the argument --gui by another internal argument --gui_launched.
        return False
    else:
        return True

def AutoGooey(fn):  # pragma: no cover
    """Automatically show a Gooey GUI if --gui is passed as the first argument, else it will just run the function as normal"""
    if check_gui_arg():
        return gooey.Gooey(fn)
    else:
        return fn


#***********************************
#                       MAIN
#***********************************


@AutoGooey
def main(argv=None, return_report=False):
    if argv is None: # if argv is empty, fetch from the commandline
        argv = sys.argv[1:]
    elif isinstance(argv, _str): # else if argv is supplied but it's a simple string, we need to parse it to a list of arguments before handing to argparse or any other argument parser
        argv = shlex.split(argv) # Parse string just like argv using shlex

    #==== COMMANDLINE PARSER ====

    #== Commandline description
    desc = '''CSG Dicoms Anonymizer v%s
Description: Anonymize dicoms and demographics using undecryptable hashs.

Note: use --cmd to avoid launching the graphical interface and use as a commandline tool.
    ''' % __version__
    ep = ''' '''

    #== Commandline arguments
    #-- Constructing the parser
    # Use GooeyParser if we want the GUI because it will provide better widgets
    if (not '--cmd' in argv and not '--ignore-gooey' in argv and not '--help' in argv and not '-h' in argv):  # pragma: no cover
        # Initialize the Gooey parser
        main_parser = gooey.GooeyParser(add_help=True, description=desc, epilog=ep, formatter_class=argparse.RawTextHelpFormatter)
        # Define Gooey widget types explicitly (because type auto-detection doesn't work quite well)
        widget_dir = {"widget": "DirChooser"}
        widget_filesave = {"widget": "FileSaver"}
        widget_file = {"widget": "FileChooser"}
        widget_text = {"widget": "TextField"}
    else: # Else in command-line usage, use the standard argparse
        # Delete the special argument to avoid unrecognized argument error in argparse
        if len(argv) > 0 and '--ignore-gooey' in argv: argv.remove('--ignore-gooey') # this argument is automatically fed by Gooey when the user clicks on Start
        if len(argv) > 0 and '--cmd' in argv: argv.remove('--cmd')
        # Initialize the normal argparse parser
        main_parser = argparse.ArgumentParser(add_help=True, description=desc, epilog=ep, formatter_class=argparse.RawTextHelpFormatter)
        # Define dummy dict to keep compatibile with command-line usage
        widget_dir = {}
        widget_filesave = {}
        widget_file = {}
        widget_text = {}

    # Required arguments
    main_parser.add_argument('-i', '--input', metavar='/some/path', type=str, required=True,
                        help='Path to the dicom root folder (dicoms will be replaced and non-dicom files deleted).', **widget_dir)
    main_parser.add_argument('-d', '--demographics', type=str, required=True, default='db.csv', metavar='db.csv',
                        help='Demographics file to anonymize.', **widget_file)

    # Optional output/copy mode
    main_parser.add_argument('--resume', action='store_true', required=False, default=False,
                        help='Resume with a previously generated anonymization map (need to provide dicom_names.csv and idtoname.csv). Can be used to generate an updated anonymized demographic file. Can also be used to resume from a bug leaving you with a partial dicoms anonymization.')
    main_parser.add_argument('--distance', type=float, required=False, default=0.2, metavar=0.2,
                        help='Distance threshold for the jaccard distance (words and letters) for similar names. normalized distance threshold to match similar names. 0.0 is no difference, 1.0 is everything different.', **widget_file)
    main_parser.add_argument('--anon_prefix', type=str, required=False, default='subj_', metavar='subj_',
                        help='Prefix of the generated anonymized ids.', **widget_text)
    main_parser.add_argument('--anon_salt', type=str, required=False, default=None, metavar='some random string',
                        help='Set this to a string of your choice to generate unique hashes, this adds protection against decrypting the ids (but keep the same if you want to be able to update anonymized data).', **widget_text)
    main_parser.add_argument('--anon_hash_algo', type=str, required=False, default='md5', metavar='md5',
                        help='Hash algorithm to generate the anonymized ids.', **widget_text)
    main_parser.add_argument('--anon_length', type=int, required=False, default=8, metavar='8',
                        help='Length of the id. Set to None to disable. Shortening might be an added security if anon_permanent_ids == True, but raises the risks of collisions (two names having same id). You can try, there will be an exception anyway if there is a collision.')
    main_parser.add_argument('--anon_irreversible_ids', action='store_true', required=False, default=False,
                        help='This adds a layer of security and also generates nice ids from 1 to number of subjects. Do not enable this if you want the anonymized id to always be the same (useful if you want to add new subjects and keep the same ids for old ones, eg, if you have multiple datasets with overlapping subjects but with different ones as well), but at the expense of security (the anonymized id can potentially be decrypted).')
    main_parser.add_argument('--dropcols', type=str, required=False, default='report_path;alt_names', metavar='item1;item2;...',  # don't forget to NOT use the | symbol, else argparse will chocke
                        help='Columns to drop from demographics file (eg, fields that can potentially contain patient name).', **widget_text)
    main_parser.add_argument('-l', '--log', metavar='/some/folder/filename.log', type=str, required=False,
                        help='Path to the log file. (Output will be piped to both the stdout and the log file)', **widget_filesave)
    main_parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False,
                        help='Verbose mode (show more output).')
    main_parser.add_argument('--silent', action='store_true', required=False, default=False,
                        help='No console output (but if --log specified, the log will still be saved in the specified file).')


    #== Parsing the arguments
    args = main_parser.parse_args(argv) # Storing all arguments to args
    
    #-- Set variables from arguments
    inputpath = get_fullpath(args.input)
    demo_csv = get_fullpath(args.demographics)
    rootpath = inputpath
    resume = args.resume
    dist_threshold = args.distance
    anon_prefix = args.anon_prefix
    anon_salt = args.anon_salt
    anon_hash_algo = args.anon_hash_algo
    anon_length = args.anon_length
    anon_permanent_ids = not args.anon_irreversible_ids
    demo_cols_drop = args.dropcols.split(';') if args.dropcols else ['report_path', 'alt_names']
    verbose = args.verbose
    silent = args.silent

    # -- Sanity checks
    if os.path.isfile(inputpath): # if inputpath is a single file (instead of a folder), then define the rootpath as the parent directory (for correct relative path generation, else it will also truncate the filename!)
        rootpath = os.path.dirname(inputpath)

    # Strip trailing slashes to ensure we correctly format paths afterward
    if rootpath:
        rootpath = rootpath.rstrip('/\\')

    if not os.path.isdir(rootpath):
        raise NameError('Specified input dicoms path does not exist. Please check the specified path')
    
    if not os.path.isfile(demo_csv) or not os.path.exists(demo_csv):
        raise NameError('Specified demographics file path does not exist. Please check the specified path')

    # -- Configure the log file if enabled (ptee.write() will write to both stdout/console and to the log file)
    if args.log:
        ptee = Tee(args.log, 'a', nostdout=silent)
        #sys.stdout = Tee(args.log, 'a')
        sys.stderr = Tee(args.log, 'a', nostdout=silent)
    else:
        ptee = Tee(nostdout=silent)

    #### Main program


    if not resume:
        # ----------------------------------------
        # # Part 1
        # ## Extract dicom names


        # In[ ]:

        # -- Get the list of dicoms (they must all be at the first level, one folder per subject)

        # Get unzipped dicom folders list
        print('Constructing list of DICOM subjects through folders names, please wait...')
        subjects_list = get_list_of_folders(rootpath)

        # Extract name from first readable dicom
        if verbose:
            print('Found subjects dicom folders: %s' % ', '.join(subjects_list))

        dcm_subj_list, folder_to_name = get_dcm_names_from_dir(rootpath, verbose=verbose)
        print('Total dicom subjects: %i. Detailed list: %s' % (len(dcm_subj_list), ', '.join(dcm_subj_list)))


        # In[ ]:

        # -- Extracting subjects names from zip files
        # NOTE: anonymization is NOT supported on zip files, only on unzipped dicom folders!

        # Extract list of zip files
        subjects_zip_list = get_list_of_zip(rootpath)
        print(subjects_zip_list)

        dcm_subj_list, folder_to_name = get_dcm_names_from_zip(rootpath, dcm_subj_list, folder_to_name, verbose=verbose)
        print('Total dicom subjects: %i. Detailed list: %s' % (len(dcm_subj_list), ', '.join(dcm_subj_list)))


        # In[ ]:

        folder_to_name


        # In[ ]:

        # Save all extracted fields to a csv file!
        output_file = 'dicom_names.csv'
        save_dict_as_csv([{'name': name, 'path': path} for name, path in zip(dcm_subj_list, subjects_list + subjects_zip_list)], output_file, csv_order_by='name', verbose=True)
        print('Dicom patients names saved to csv file: %s' % output_file)


        # ---------------------------------
        # # Part 2: Generate anonymization mapping
        # ## Anonymization initialization (generate anonymized ids)


        # In[ ]:
        with open("dicom_names.csv") as f:
            dcm_subj_list = [row['name'] for row in csv.DictReader(f, delimiter=';')]
        print(dcm_subj_list)


        # In[ ]:

        # Generate disambiguated list of dicom names
        # Disambiguate dicom names
        cd_unique = disambiguate_names(dcm_subj_list, verbose=True)
        # Extract list of unique dicom names
        dcm_unique = set([c['name'] for c in cd_unique])

        print(dcm_unique)


        # In[ ]:

        # Unique id extractor from name, insensitive to accentuated charaters nor non-alphabetical characters nor firstname/lastname position switching
        from csg_fileutil_libs.aux_funcs import sort_list_a_given_list_b, replace_buggy_accents, _unidecode

        def clean_name(name):
            '''Clean name from accents and non-alphabetical characters'''
            return re.sub(r'\W', r'', _unidecode(replace_buggy_accents(name.decode('utf8'), 'utf8')).lower())
        def extract_ordered_letters(name):
            '''Order letters composing a name to alphabetical order'''
            alphabet = list('abcdefghijklmnopqrstuvwxyz1234567890-')
            return ''.join(sort_list_a_given_list_b(list(clean_name(name)), alphabet))
        def get_hash(string, algo=None):
            if algo is None or algo == 'md5':
                return hashlib.md5(string).hexdigest()
            elif algo == 'sha1':
                return hashlib.md5(string).hexdigest()
            else:
                raise NameError('Hash algorithm not recognized: %s' % algo)
        def get_ordered_hash(hash_func, string, salt=None, algo=None):
            '''Get a unique hash insensitive to accents, non-alphabetical characters nor words position switching'''
            return hash_func(extract_ordered_letters(string+(salt if salt else '')), algo=algo)
        def get_duplicates(d):
            seen = set()
            for k, v in d.items():
                if v in seen:
                    yield k, v
                else:
                    seen.add(v)

        # Unit test
        name1 = 'rajaé chatila'
        name2 = 'chatila  rajaé|'
        assert(get_ordered_hash(get_hash, name1) == get_ordered_hash(get_hash, name2))


        # In[ ]:

        ##### Generate anonymization scheme from dicoms patients names #####

        # Generate unique hashes from each dicom's patient name
        # Generate an anonymized id resilient to spaces and non letters characters and words switching
        # to do that, we take the name, and reorder all letters (and remove any non-letter symbol) by alphabetical order, which gives us simply the ordered sequence of letters composing each name
        anon_hashes = {name: get_ordered_hash(get_hash, name, anon_salt, algo=anon_hash_algo) for name in dcm_unique}
        # Shortening is an added security, so that if someone tries to bruteforce, there will be missing info to reconstitute the original name that gave this hash (because we are missing parts of the hash, so lots of dissimilar names will have the same shortened hash)
        if anon_length:
            for name, h in anon_hashes.items():
                anon_hashes[name] = anon_hashes[name][:anon_length]
        # There can be collisions in hashes, then check that there is none
        anon_dups = dict(get_duplicates(anon_hashes))
        if anon_dups:
            anon_dups_print = {h: [name for name, h2 in anon_hashes.items() if h2 == h] for h in anon_dups.values()}
            raise ValueError('Two names have the same id! Please use another hashing algorithm or raise hash length or another salt or turn off anon_permanent_ids. Here is the list of names with same ids: %s' % anon_dups_print)
        del anon_dups

        # Generate the final id (second step)
        if anon_permanent_ids:
            # generate a straightforward id from a shortened hash
            names_and_ids = anon_hashes
        else:
            # generate a unique id based on order (simply the order number when ordered by the hash - since we use the "ordered hash", we get the same properties: the same set of patients names will always generate the same order, and the order cannot be traced back, since it depends both on the hash AND the exact set of patients names to get the exact same order)
            # the big advantage of this approach is that it is nearly impossible to decrypt the original name, since the id gives strictly no information at all
            # the disadvantage is that it is dependent on the subjects names list, so if you add a subject, nearly all ids will change
            names_and_ids = {name: str(id+1).zfill(anon_length) for id, name in enumerate(sorted(anon_hashes, key=anon_hashes.get))}
        # Prepend prefix and save the anonymized ids
        anon_ids = {("%s%s" % (anon_prefix, id)): name for name, id in names_and_ids.items()}
        anon_ids


        # In[ ]:

        # Write anonymization csv
        ids_to_name = [{'id': pts_id, 'name': anon_ids[pts_id]} for pts_id in sorted(anon_ids)]
        save_dict_as_csv(ids_to_name, 'idtoname.csv', fields_order=['id', 'name'], csv_order_by='id', verbose=False)
        print('Conversion list (id -> name) saved to idtoname.csv.')


    # ------------------------------------
    # # Part 3: applying anonymization
    # TODO: remake to adapt the dcm_names to the ones in the specified rootpath, not the ones in dicom_names.csv!
    # ## Merging dicom names and demographics csv names


    # In[ ]:
    with open(demo_csv) as f:
        cf = list(csv.DictReader(f, delimiter=';'))

    with open("dicom_names.csv") as f:
        dcm_subj_list = [row['name'] for row in csv.DictReader(f, delimiter=';')]


    # In[ ]:

    # Generate disambiguated list of dicom names
    # Disambiguate dicom names
    cd_unique = disambiguate_names(dcm_subj_list, verbose=verbose)
    # Extract list of unique dicom names
    dcm_unique = set([c['name'] for c in cd_unique])

    print(dcm_unique)

    # Generate dicom name to unique name mapping
    dcmname_to_uniquename = {(c['alt_names'] if 'alt_names' in c else c['name']): c['name'] for c in cd_unique}

    print(dcm_unique)
    print(dcmname_to_uniquename)


    # In[ ]:

    # Disambiguate and clean up csv names
    # Cleanup names
    for i in range(len(cf)):
        cf[i]['name'] = cleanup_name(cf[i]['name'])

    # Disambiguate (ie, same name with typos or inversed firstname/lastname)
    cf = disambiguate_names(cf, dist_threshold=dist_threshold, verbose=verbose)
    # Print list of disambiguated names
    [{c['name']: c['alt_names']} for c in cf if 'alt_names' in c and c['alt_names']]


    # In[ ]:

    # Computing distance matrix (ie, finding similar names between dicoms and demographics csv)
    print('Computing distance matrix (finding similar names) between dicoms and demographics, please wait...')
    dist_matches = {}
    name_to_anon_ids = {v: k for k, v in anon_ids.items()}
    for subj in _tqdm(dcm_unique, desc='distmat', unit='subj', file=sys.stdout):
        found = False
        for c in cf:
            if distance.nlevenshtein(subj, c['name'], method=1) <= dist_threshold or distance_jaccard_words_split(subj, c['name'], partial=True, norm=True, dist=dist_threshold) <= dist_threshold: # use shortest distance with normalized levenshtein
                if subj not in dist_matches:
                    dist_matches[subj] = []
                dist_matches[subj].append(c['name'])
                found = True
        if not found:
            dist_matches[subj] = None

    # Remove duplicate values (ie, csv names)
    dist_matches = {k: (list(set(v)) if v else v) for k, v in dist_matches.items()}
    # Find missing subjects (ie, dicom name present but missing in csv database)
    missing_subj = {k: v for k, v in dist_matches.items() if not v}
    # Print results
    if not missing_subj:
        print('No missing subject, congratulations!')
    else:
        print('Missing subjects from csv database (saved in missing_demo.csv): %i, names: %s' % (len(missing_subj), ', '.join(sorted(missing_subj.keys()))))
        save_dict_as_csv([{'name': msubj, 'id': name_to_anon_ids[msubj]} for msubj in missing_subj.keys()], 'missing_demo.csv', fields_order=['name', 'id'], csv_order_by='name', verbose=False)
        save_dict_as_csv([{'id': name_to_anon_ids[msubj]} for msubj in missing_subj.keys()], 'missing_demo_anonymized.csv', fields_order=['id'], csv_order_by='id', verbose=False)

    print('\nList of all matches (dicom : csv):')
    print(dist_matches)


    # In[ ]:

    # Compute the csv name to dicom unique name mapping
    csvname_to_uniquename = {}
    for uniquename in dist_matches.keys():
        if dist_matches[uniquename]:
            for csv_name in dist_matches[uniquename]:
                csvname_to_uniquename[csv_name] = uniquename


    # In[ ]:

    def flatten_gen(L):
        for item in L:
            if isinstance(item, list):
                for i in flatten(item):
                    yield i
            else:
                yield item

    def flatten(L):
        return list(flatten_gen(L))

    def get_unique_names(L):
        return filter(None, set(L))


    # In[ ]:

    demo_short_csv = 'demographics_shortened.csv'
    # Get unique demo csv names (which matched with dicom in the distance matrix)
    demo_names = get_unique_names(flatten(dist_matches.values()))
    # Shorten demographics to only names present in dicoms
    cf_short = [c for c in cf if c['name'] in demo_names]
    # Save shortened demographics
    save_dict_as_csv(cf_short, demo_short_csv, fields_order=['name'], csv_order_by='name', verbose=False)
    print('Shortened demographics (to only the dicoms available) were saved to %s.' % demo_short_csv)


    # In[ ]:

    cf_short


    # -------------------------------
    # ## Anonymizing demographics csv

    # In[ ]:
    with open('idtoname.csv', mode='r') as f:
        reader = csv.DictReader(f, delimiter=';')
        anon_ids = {row['id']: row['name'] for row in reader}
    anon_ids


    # In[ ]:

    name_to_anon_ids = {v: k for k, v in anon_ids.items()}
    name_to_anon_ids


    # In[ ]:

    # Anonymize demographics csv

    demo_anon_csv = 'demographics_anonymized.csv'

    # Make a new dataframe from shortened demographics csv
    cf_anon = list(cf_short)  # copy
    # Anonymize names by using the csv name -> unique name -> anonymized id mapping
    # TODO: might be better to use pandas join?
    for rowid in range(len(cf_anon)):
        cf_anon[rowid]['name'] = name_to_anon_ids[csvname_to_uniquename[cf_anon[rowid]['name']]]
    # Drop columns that we cannot anonymize but might give away subjects infos (like report_path and alt_names)
    if demo_cols_drop:
        for rowid in range(len(cf_anon)):
            for col in demo_cols_drop:
                if col in cf_anon[rowid]:
                    del cf_anon[rowid][col]
    # Save anonymized demographics
    save_dict_as_csv(cf_anon, demo_anon_csv, fields_order=['name'], csv_order_by='name', verbose=False)
    print('Anonymized demographics successfully saved to %s.' % demo_anon_csv)


    # -----------------
    # ## Anonymizing dicoms

    # In[ ]:

    # Unzip all into folders, because we cannot anonymize zip files (ie, can't modify files inside a zip)
    # NOTE: ZIP FILES WILL BE DELETED!
    def unzip(zipfilepath, outputpath):
        zip_ref = zipfile.ZipFile(zipfilepath, 'r')
        zip_ref.extractall(outputpath)
        zip_ref.close()

    # First we need to delete all .DS_Store files (else we get permission denied IOError)
    count_dstore = 0
    for dirpath, filename in _tqdm(recwalk(rootpath, topdown=False, folders=True), unit='files', desc='DSTOREDEL', file=sys.stdout):
        if filename.lower() == '.ds_store':
            fullfilepath = os.path.join(dirpath, filename)
            os.remove(fullfilepath)
            count_dstore += 1
    print('Total .DS_Store files deleted: %i.' % count_dstore)

    # Unzip files and delete zip
    subjects_list_zip = get_list_of_zip(rootpath)
    count_zip = 0
    for z in _tqdm(subjects_list_zip, unit='files', desc='UNZIP', file=sys.stdout):
        zipfilepath = os.path.join(rootpath, z)
        if verbose:
            print('- Unzipping file: %s' % z)
        # Unzip file into the same root directory as other dicom folders
        unzip(zipfilepath, rootpath)
        # Delete the zip (since we cannot anonymize it)
        os.remove(zipfilepath)
        count_zip += 1
    print('Total zipfiles unzipped: %i.' % count_zip)


    # In[ ]:

    # -- Anonymization of dicom files
    # Note: this will anonymize only the dicoms fields (name of patient, whatever field it is found in).
    # If there are any other file containing the patient's name (such as .txt, .csv, .xls, etc), the files might be deleted if you want (add the extension in the list) or they will stay.
    #from dicom.filebase import DicomFileLike  # fix for IOError access denied, see https://github.com/darcymason/pydicom/issues/69
    def find_hidden_name_fields(dcmdata, dcm_pts_names, hidden_name_fields=None):
        '''From a pydicom object, return all fields where one of the dcm_pts_name (a list) is present.
        This ease the detection of additional fields where patient name was stored.'''
        if hidden_name_fields is None:
            hidden_name_fields = set()
        # Convert name to regex friendly (because dicoms often replace spaces by ^)
        dcm_pts_names = [pts_name.replace(' ', '[\W]+') for pts_name in dcm_pts_names]
        # Walk through each dicom field
        for dcmfield in dcmdata.keys(): # different from dir()?
            if dcmdata[dcmfield]:
                try:
                    #dcmfieldval = dcmdata.data_element(dcmfield).value
                    dcmfieldval = dcmdata[dcmfield].value
                    check = False
                    if isinstance(dcmfieldval, list):
                        dcmfieldval_lower = [s.lower() if isinstance(s, str) else s for s in dcmfieldval]
                        check = any(pts_name in dcmfieldval_lower for pts_name in dcm_pts_names)
                    elif isinstance(dcmfieldval, (int, float)):
                        check = False
                    else:
                        check = any(re.search(pts_name, dcmfieldval.lower()) for pts_name in dcm_pts_names)
                    if check:
                        hidden_name_fields.add(dcmfield)
                except AttributeError:
                    print('Error with field: %s' % str(dcmfield))
                    raise
        return hidden_name_fields


    reports_delete = True  # delete pdf/doc/docx/txt files automatically?
    skip_already_processed = True
    remove_private_tags = False
    fields_to_del = ['PatientAddress', 'PatientBirthTime', 'PatientTelephoneNumbers', 'OtherPatientNames']
    print('-- Anonymization started, please wait, this might take a while (also make sure you unzipped all dicoms into folders)...')
    print('Note: if you get an IOError permission denied error, make sure you close any file explorer or application using any of the subjects folder (including Windows Explorer, else folders cannot be renamed).')
    print('Note2: JPEG2000 compressed dicom files are unsupported, please uncompress them beforehand (eg, using dcmdjpeg).')
    print('Note3: in case of an Access Error, you can continue the anonymization, it will restart from the start but it will skip already processed dicom files.')
    count_anon = 0
    count_files = 0
    count_delete = 0
    count_files_skipped = 0
    # Init path and 1st level folders list
    uni_rootpath = unicode(rootpath, 'latin1')  # convert rootpath to unicode before walking with os.listdir and recwalk, so we get back unicode strings too (else we won't be able to enter folders with accentuated characters)
    subjects_list = get_list_of_folders(uni_rootpath)
    # Precompute total number of files (for progressbar)
    print('Precomputing total number of files, please wait...')
    for subject in _tqdm(subjects_list, unit='folders', desc='PRECOMP', file=sys.stdout):
        fullpath = os.path.join(uni_rootpath, subject)
        for dirpath, filename in recwalk(fullpath, topdown=False, folders=True):
            count_files += 1
    # Get folder_to_name mapping
    _, folder_to_name = get_dcm_names_from_dir(uni_rootpath)
    _, folder_to_name = get_dcm_names_from_zip(uni_rootpath, folder_to_name=folder_to_name)
    # Loop through each subject root directory to rewrite dicoms
    print('Launching anonymization of dicoms fields, please wait...')
    tbar = _tqdm(total=count_files, unit='files', desc='ANON', file=sys.stdout)
    hfields = set()
    for folder in subjects_list:
        # Already processed folder and there are several sessions, extract the id from folder name
        subject = folder
        try:
            subject = re.match('(^%s.+)_s\d+$' % anon_prefix, subject).group(1)
        except AttributeError:
            pass
        # Already processed folder, then retrieve back the patient's name from the anonymized id
        if subject in anon_ids:
            pts_name = anon_ids[subject]
            if skip_already_processed:
                fullpath = os.path.join(uni_rootpath, folder)
                c = 0
                for _, _ in recwalk(fullpath, topdown=False, folders=True):
                    c += 1
                tbar.update(c)
                count_files_skipped += c
                continue
        # Partially processed folder, get the original name and continue
        elif folder_to_name[subject] in anon_ids:
            pts_name = anon_ids[folder_to_name[subject]]
        else:
            pts_name = folder_to_name[subject]
        # Special case: no dicom with a patient name can be found inside the folder (might be nifti files instead?), so we just skip
        if pts_name is None:
            continue
        # Fetch the anonymized id from folder name (because we already looked inside to get the first dicom's patientname)
        anon_id = name_to_anon_ids[dcmname_to_uniquename[pts_name]]
        if verbose:
            print('- Processing subject %s -> %s in folder %s' % (pts_name, anon_id, folder))
        #fullpath = unicode(os.path.join(rootpath, folder), 'latin1')
        fullpath = os.path.join(uni_rootpath, folder)  # no need to use unicode(str, 'latin1') here because rootpath and folder were both converted to unicode before
        # Loop through each subfiles and subfolders for this subject (we assume all dicoms are for one subject, so we rename them all to this subject)
        for dirpath, filename in recwalk(fullpath, topdown=False, folders=True):
            fullfilepath = os.path.join(dirpath, filename)
            # Report file: delete if option enabled
            if reports_delete and filename.endswith( ('pdf', 'doc', 'docx', 'txt', 'csv', 'xls', 'xlsx') ):
                os.remove(fullfilepath)
                count_delete += 1
                continue
            elif os.path.isdir(fullfilepath):  # else we get an IOError...
                continue
            else:
                # Dicom file: change PatientName field
                try:
                    #TODO: autodetect if name is in filename and change!
                    #print('* Try to read fields from dicom file: %s' % os.path.join(dirpath, filename))
                    # Read dicom's file data
                    #os_id = os.open(str(fullfilepath), os.O_BINARY | os.O_RDONLY)
                    #fd = os.fdopen(os_id)
                    #dcmdata = dicom.read_file(DicomFileLike(fd), stop_before_pixels=False)
                    dcmdata = dicom.read_file(fullfilepath, stop_before_pixels=False)  # need to read the full dicom here since we will modify it, so stop_before_pixels must be False
                    # Store current name (to check at the end if we correctly cleaned up the name)
                    try:
                        dcm_pts_name = _unidecode(dcmdata.PatientName.decode('latin1').replace('^', ' ')).lower().strip()
                        # Already anonymized dicom? Get the original patient's name from the anonymized id
                        if dcm_pts_name in anon_ids:
                            dcm_pts_name = anon_ids[dcm_pts_name]
                            if skip_already_processed:
                                tbar.update()
                                continue
                    except AttributeError as exc:
                        if filename.upper() == 'DCMDIR' or filename.upper() == 'DICOMDIR':
                            os.remove(fullfilepath)  # DICOMDIR files are useless, they are only descriptive files for CD/DVD of dicoms
                            # DOES NOT WORK: pydicom can read and edit dicomdir files but cannot save them yet!
                            #dcmdata = dicom.read_dicomdir(r'dicoms\ANTOINE_el\EPI_T1\DICOMDIR')
                            #for record in dcmdata.patient_records:
                                #record.PatientName = anon_id
                            continue
                        else:
                            raise
                    # Anonymize
                    dcmdata.PatientName = anon_id
                    dcmdata.PatientID = anon_id
                    if [0x33,0x1013] in dcmdata:  # custom patientname field...
                        dcmdata[0x33,0x1013].value = anon_id
                    # Delete private fields
                    for field in fields_to_del:
                        if field in dcmdata:
                            if isinstance(field, str):
                                del dcmdata[dcmdata.data_element(field).tag]
                            else:
                                del dcmdata[field]
                    if remove_private_tags:
                        dcmdata.remove_private_tags()
                    # Try to anonymize hidden name fields
                    hfields = find_hidden_name_fields(dcmdata, [dcm_pts_name, pts_name], hfields)
                    for dcmfield in hfields:
                        if dcmfield in dcmdata:
                            dcmdata[dcmfield].value = re.sub(dcm_pts_name.replace(' ', '[\W]+'), anon_id, dcmdata[dcmfield].value, flags=re.I)
                            dcmdata[dcmfield].value = re.sub(pts_name.replace(' ', '[\W]+'), anon_id, dcmdata[dcmfield].value, flags=re.I)
                    # Last check just in case we could not remove the name everywhere!
                    dcm_data_str = _unidecode(str(dcmdata).decode('latin1').replace('^', ' ')).lower().strip()  # read all fields at once
                    if dcm_pts_name in str(dcm_data_str) or pts_name in str(dcm_data_str):  # if patient's name is still in the file, that's bad!
                        print('Hidden name fields found: %s' % hfields)
                        print('names: %s - %s' %(dcm_pts_name, pts_name))  # debugline
                        print(str(dcm_data_str))
                        raise ValueError('Error: could not remove name totally (there must be an additional non-standard PatientName field) from file: %s' % fullfilepath)
                    # Save anonymized dicom file
                    dcmdata.save_as(fullfilepath)
                    # Close the dicom file
                    #os.close(os_id)
                    del dcmdata
                    count_anon += 1
                except (InvalidDicomError) as exc:
                    pass
                except AttributeError as exc:
                    print(fullfilepath)
                    raise
            tbar.update()  # update progressbar
    tbar.close()

    print('Hidden name fields found (and automagically anonymized): %s' % hfields)
    print('Total dicom anonymized: %i over %i total. Total dicom files skipped: %i. Total reports/non-dicom files deleted: %i.' % (count_anon, count_files, count_files_skipped, count_delete))


    # In[ ]:

    # Rename files if filename include a patient's name

    # Compile regex to find any patient name (of any patient!) in a string. Non-alphabetical characters are ignored.
    filename_patterns = re.compile('(' + '|'.join(re.sub('[^a-zA-Z]+', '[^a-zA-Z]*', s) for s in dcm_unique) + ')', flags=re.I)

    uni_rootpath = unicode(rootpath, 'latin1')  # convert rootpath to unicode before walking with os.listdir and recwalk, so we get back unicode strings too (else we won't be able to enter folders with accentuated characters)
    subjects_list = get_list_of_folders(uni_rootpath)

    count_files = 0
    # Precompute total number of files (for progressbar)
    print('Precomputing total number of files, please wait...')
    for subject in _tqdm(subjects_list, unit='folders', desc='PRECOMP', file=sys.stdout):
        fullpath = os.path.join(uni_rootpath, subject)
        for dirpath, filename in recwalk(fullpath, topdown=False, folders=True):
            count_files += 1

    # Rename files if they have a patient's name
    count_moved = 0
    tbar = _tqdm(total=count_files, unit='files', desc='ANONFN', file=sys.stdout)
    print('Anonymizing of file/folder names, please wait...')
    for folder in subjects_list:  # do not rename the top directories, this will be done separately
        if verbose:
            print('- Processing top folder %s' % (folder))
        fullpath = os.path.join(uni_rootpath, folder)
        for dirpath, filename in recwalk(fullpath, topdown=False, folders=True):
            # Find any name (of any patient) in the filename
            # TODO: construct re all permutations of all names, and re.compile, it will be fast
            # TODO: try to do levenshtein distance on names? (but just with current patient name, else it will take too much time with all patients...) it will considerably slow down the anonymization... Is there a faster way?
            matchs = filename_patterns.finditer(filename)
            # If found, we find the anonymized id for each match to replace
            to_replace = []
            for m in matchs:
                # Clean up the name
                pts_name_in_filename = re.sub('[^a-zA-Z]+', ' ', m.group(1).lower())
                # Find the closest unique name
                dst_mat = dist_matrix([pts_name_in_filename], dcm_unique)
                # Get the anonymized id from unique name
                if dst_mat[pts_name_in_filename]:
                    anon_id = name_to_anon_ids[dst_mat[pts_name_in_filename][0]]
                else:  # could not find an id, just anonymize with a random name
                    anon_id = 'anon'
                # Add slide index and anonymized id to replace all at once later
                to_replace.append( ( anon_id, slice(m.start(1), m.end(1)) ) )
            # Replace all matchs at once
            if to_replace:
                # Can't modify strings, need to convert to a list
                filename_anon = list(filename)
                # For each match, replace with anonymized id
                for anon_id, slidx in to_replace[::-1]:  # reverse list because else the subsequent items won't be aligned anymore when we will replace the first items in the list
                    filename_anon[slidx] = anon_id
                # Convert back to a string
                filename_anon = ''.join(filename_anon)
                # Rename the file/folder
                shutil.move(os.path.join(dirpath, filename), os.path.join(dirpath, filename_anon))
                count_moved += 1
            tbar.update()
    tbar.close()
    print('Total dicom files/folders moved: %i over %i total.' % (count_moved, count_files))


    # In[ ]:

    # Rename folders if enabled
    rename_folders = True  # rename root folder # TODO: rename any subfolder with patient's name
    delete_empty = True  # delete empty folders (or not containing any DICOM, such as nifti folders)?

    count_folder = 0
    count_skipped = 0
    count_empty = 0
    uni_rootpath = unicode(rootpath, 'latin1')  # convert rootpath to unicode before walking with os.listdir and recwalk, so we get back unicode strings too (else we won't be able to enter folders with accentuated characters)
    if rename_folders:
        print('Launching anonymization of dicom folders, please wait...')
        # Get list of folders
        subjects_list = get_list_of_folders(uni_rootpath)
        for subject in _tqdm(subjects_list, unit='folder', desc='RENAME', file=sys.stdout):
            # Already anonymized folder, just skip
            if subject in anon_ids or re.match('(^%s.+)_s\d+$' % anon_prefix, subject):
                count_skipped += 1
                continue
            pts_name = folder_to_name[subject]
            # Already partially anonymized, we get the original name from the anonymization mapping
            if pts_name in anon_ids:
                pts_name = anon_ids[pts_name]
            # Special case: no dicom with a patient name can be found inside the folder (might be nifti files instead?), so we just remove this folder!
            if pts_name is None and delete_empty:
                filepath = os.path.join(uni_rootpath, subject)
                if os.path.isdir(filepath):
                    shutil.rmtree(filepath)
                else:
                    os.remove(filepath)
                count_empty += 1
                continue
            # Fetch the anonymized id from folder name (because we already looked inside to get the first dicom's patientname)
            anon_id = name_to_anon_ids[dcmname_to_uniquename[pts_name]]
            if verbose:
                print('- Processing subject %s -> %s in folder %s' % (pts_name, anon_id, subject))
            fullpath = os.path.join(uni_rootpath, subject)
            # Rename subject directory
            new_folder_name = os.path.join(uni_rootpath, anon_id)
            if not os.path.exists(new_folder_name):
                os.rename(fullpath, new_folder_name)
            else:
                # if new folder already exists, find a new name (append "_sx" where x is a number)
                for i in range(2, 1000):
                    alt_folder_name = "%s_s%i" % (new_folder_name, i)
                    if not os.path.exists(alt_folder_name):
                        os.rename(fullpath, alt_folder_name)
                        break
            count_folder += 1

    print('Total dicom folders renamed (anonymized): %i over %i total. Skipped: %i. Empty folders (or containing non-dicom files) and thus deleted: %i.' % (count_folder, len(subjects_list), count_skipped, count_empty))


    # In[ ]:

    # Shorten (again) anonymized demographics to only the subjects we have dicom folders for
    demo_anon_csv = 'demographics_anonymized.csv'
    with open(demo_anon_csv) as f:
        cf_anon = list(csv.DictReader(f, delimiter=';'))
    # Get list of anonymized dicom names
    dcm_ids, _ = get_dcm_names_from_dir(rootpath)
    # Shorten anonymized demographics to only the ids present in dicoms
    for rowid in range(len(cf_anon)):
        if cf_anon[rowid]['name'] not in dcm_ids:
            del cf_anon[rowid]
            for col in demo_cols_drop:
                if col in cf_anon[rowid]:
                    del cf_anon[rowid][col]
    # Save anonymized demographics
    demo_anon_short_csv = 'demographics_anonymized_shortened.csv'
    save_dict_as_csv(cf_anon, demo_anon_short_csv, fields_order=['name'], csv_order_by='name', verbose=False)
    print('Shortened anonymized demographics (to only the dicoms available) were saved to %s.' % demo_anon_short_csv)

    # TODO: PANDAS VERSION: to delete if the above works (no pandas = smaller standalone package)
    # Shorten (again) anonymized demographics to only the subjects we have dicom folders for
    #demo_anon_csv = 'demographics_anonymized.csv'
    #cf_anon = pd.read_csv(demo_anon_csv, sep=';').fillna('')
    # Get list of anonymized dicom names
    #dcm_ids, _ = get_dcm_names_from_dir(rootpath)
    # Shorten anonymized demographics to only the ids present in dicoms
    #cf_anon = cf_anon[cf_anon['name'].isin(dcm_ids)]
    # Save shortened anonymized demographics
    #cf_anon.to_csv(demo_anon_csv, sep=';', na_rep='NA', index=False)
    #print('Shortened anonymized demographics (to only the dicoms available) were saved to %s.' % demo_anon_csv)


    # In[ ]:

    # TODO: add new ids (new dicom folders) even if not in demographics
    
    return 0


# Calling main function if the script is directly called (not imported as a library in another program)
if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
