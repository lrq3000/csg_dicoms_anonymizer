#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Auxiliary functions library for reports extractor and dicom anonymization
# Copyright (C) 2017 Larroque Stephen
# Licensed under MIT License.
#

from __future__ import absolute_import

import csv
import os
import re
from .distance import distance

try:
    from scandir import walk # use the faster scandir module if available (Python >= 3.5), see https://github.com/benhoyt/scandir
except ImportError as exc:
    from os import walk # else, default to os.walk()

try:
    # to convert unicode accentuated strings to ascii
    from .unidecode import unidecode
    _unidecode = unidecode
except ImportError as exc:
    # native alternative but may remove quotes and some characters (and be slower?)
    import unicodedata
    def _unidecode(s):
        return unicodedata.normalize('NFKD', s).encode('ascii', 'ignore')
    print("Notice: for reliable ascii conversion, you should pip install unidecode. Falling back to native unicodedata lib.")

try:
    from .tqdm import tqdm
    _tqdm = tqdm
except ImportError as exc:
    def _tqdm(*args, **kwargs):
        if args:
            return args[0]
        return kwargs.get('iterable', None)

def save_dict_as_csv(d, output_file, fields_order=None, csv_order_by=None, verbose=False):
    """Save a dict/list of dictionaries in a csv, with each key being a column"""
    # Define CSV fields order
    # If we were provided a fields_order list, we will show them first, else we create an empty fields_order
    if fields_order is None:
        fields_order = []
    # Get dict/list values
    if isinstance(d, dict):
        dvals = list(d.values())
    else:
        dvals = d
    # dict is empty, maybe no match was found? Then we just save an empty csv
    if not dvals:
        with open(output_file, 'wb') as f:
            f.write('')
        return
    # Then automatically add any other field (which order we don't care, they will be appended in alphabetical order)
    fields_order_check = set(fields_order)
    for missing_field in sorted(dvals[0]):
        if missing_field not in fields_order_check:
            fields_order.append(missing_field)
    if verbose:
        print('CSV fields order: '+str(fields_order))

    # Write the csv
    with open(output_file, 'wb') as f:  # Just use 'w' mode in 3.x
        w = csv.DictWriter(f, fields_order, delimiter=';')
        w.writeheader()
        # Reorder by name (or by any other column)
        if csv_order_by is not None:
            d_generator = sorted(dvals, key=lambda x: x[csv_order_by])
        else:
            d_generator = dvals
        # Walk the ordered list of dicts and write each as a row in the csv!
        for d_fields in d_generator:
            w.writerow(d_fields)
    return True


def save_df_as_csv(d, output_file, fields_order=None, csv_order_by=None, verbose=False):
    """Save a dataframe in a csv"""
    # Define CSV fields order
    # If we were provided a fields_order list, we will show them first, else we create an empty fields_order
    if fields_order is None:
        fields_order = []
    # Then automatically add any other field (which order we don't care, they will be appended in alphabetical order)
    fields_order_check = set(fields_order)
    for missing_field in sorted(d.columns):
        if missing_field not in fields_order_check:
            fields_order.append(missing_field)
    if verbose:
        print('CSV fields order: '+str(fields_order))

    # Write the csv
    d.sort_values(csv_order_by).to_csv(output_file, sep=';', index=False, columns=fields_order)
    return True


def distance_jaccard_words(seq1, seq2, partial=True, norm=False, dist=0, minlength=0):
    """Jaccard distance on two lists of words. Any permutation is tested, so the resulting distance is insensitive to words order."""
    # The goal was to have a distance on words that 1- is insensible to permutation ; 2- returns 0.2 or less if only one or two words are different, except if one of the lists has only one entry! ; 3- insensible to shortened name ; 4- allow for similar but not totally exact words.
    seq1_c = filter(None, list(seq1))
    seq2_c = filter(None, list(seq2))
    count_total = len(seq1_c) + len(seq2_c)
    count_eq = 0
    for s1 in seq1_c:
        flag_eq = False
        for skey, s2 in enumerate(seq2_c):
            if minlength and (len(s1) < minlength or len(s2) < minlength):
                continue
            if s1 == s2 or \
            (partial and (s1.startswith(s2) or s2.startswith(s1))) or \
            (dist and distance.nlevenshtein(s1, s2, method=1) <= dist):
                count_eq += 1
                del seq2_c[skey]
                break
    # Prepare the result to return
    if norm is None:
        # Just return the count of equal words
        return count_eq
    else:
        count_eq *= 2  # multiply by two because everytime we compared two items
        if norm:
            # Normalize distance
            return 1.0 - (float(count_eq) / count_total)
        else:
            # Return number of different words
            return count_total - count_eq

def distance_jaccard_words_split(s1, s2, *args, **kwargs):
    """Split sentences in words and call distance jaccard for words"""
    wordsplit_pattern = kwargs.get('wordsplit_pattern', None)
    if 'wordsplit_pattern' in kwargs:
        del kwargs['wordsplit_pattern']
    if not wordsplit_pattern:
        wordsplit_pattern = r'-+|\s+|,+|\.+|/+';

    return distance_jaccard_words(re.split(wordsplit_pattern, s1), re.split(wordsplit_pattern, s2), *args, **kwargs)

def fullpath(relpath):
    '''Relative path to absolute'''
    if (type(relpath) is object or hasattr(relpath, 'read')): # relpath is either an object or file-like, try to get its name
        relpath = relpath.name
    return os.path.abspath(os.path.expanduser(relpath))

def recwalk(inputpath, sorting=True, folders=False, topdown=True, filetype=None):
    '''Recursively walk through a folder. This provides a mean to flatten out the files restitution (necessary to show a progress bar). This is a generator.'''
    if filetype and isinstance(filetype, list):
        filetype = tuple(filetype)  # str.endswith() only accepts a tuple, not a list
    # If it's only a single file, return this single file
    if os.path.isfile(inputpath):
        abs_path = fullpath(inputpath)
        yield os.path.dirname(abs_path), os.path.basename(abs_path)
    # Else if it's a folder, walk recursively and return every files
    else:
        for dirpath, dirs, files in walk(inputpath, topdown=topdown):	
            if sorting:
                files.sort()
                dirs.sort()  # sort directories in-place for ordered recursive walking
            # return each file
            for filename in files:
                if not filetype or filename.endswith(filetype):
                    yield (dirpath, filename)  # return directory (full path) and filename
            # return each directory
            if folders:
                for folder in dirs:
                    yield (dirpath, folder)

def sort_list_a_given_list_b(list_a, list_b):
    return sorted(list_a, key=lambda x: list_b.index(x))

def replace_buggy_accents(s, encoding=None):
    """Fix weird encodings that even ftfy cannot fix"""
    # todo enhance speed? or is it the new regex on name?
    dic_replace = {
        '\xc4\x82\xc2\xa8': 'e',
        'ĂŠ': 'e',
        'Ăť': 'u',
        'â': 'a',
        'Ă´': 'o',
        'Â°': '°',
        'â': "'",
        'ĂŞ': 'e',
        'ÂŤ': '«',
        'Âť': '»',
        'Ă': 'a',
        'AŠ': 'e',
        'AŞ': 'e',
        'A¨': 'e',
        'A¨': 'e',
        'Ă': 'E',
        'â˘': '*',
        'č': 'e',
        '’': '\'',
    }
    for pat, rep in dic_replace.items():
        if encoding:
            pat = pat.decode(encoding)
            rep = rep.decode(encoding)
        s = s.replace(pat, rep)
    return s

def cleanup_name(s, encoding='latin1'):
    s = _unidecode(s.decode(encoding).replace('^', ' ')).lower().strip()
    s = re.sub('\-+', '-', re.sub('\s+', ' ', re.sub('[^a-zA-Z0-9\-]', ' ', s))).strip().replace('\r', '').replace('\n', '').replace('\t', '').replace(',', ' ').replace('  ', ' ').strip()  # clean up spaces, punctuation and double dashes in name
    return s
