csg_dicoms_anonymizer
======================

|PyPI-Status| |PyPI-Versions| |LICENCE|

This dicoms anonymizer was made to make your life easier as a clinician sharing data. It runs under Python 2.7 (it was not tested nor developped with compatibility with Python 3 in mind, although it might work with some slight changes).

This is a generic dicoms anonymizer with an emphasis on:

1. **Ease of use** : simply point to the folder containing all your dicoms (folders or zip files) and to the demographics file (just ensure you have a "name" column), and click "Start"! Note that you can also use it in commandline, allowing batch scripting.

2. **Easy to anonymize demographics** : a fuzzy matching algorithm will take care of small typos and of firstname/lastname reversal (works also if there are several first names or last names). At the end, several csv files will be generated: the anonymized demographics shortened to only the patients present in the dicoms (so that you can reuse a big demographics spreadsheet of all your subjects even to anonymize only a few subjects, the shortened anonymized demographics will only contain the pertinent subjects), the list of missing subjects that are in dicoms but not in demographics (so that your recipient knows directly what demographics were not available), and a mapping to convert from anonymized ids to original patient name (useful if your collaborator requests additional information about an anonymized patient, or just to check if there is an issue).

3. **Retaining potentially pertinent informations in dicoms** : instead of anonymizing all personal fields, this application anonymizes only the name (and strips out some other personal infos such as phone number, patient address, etc), but not all personal fields. Thus, anonymized dicoms will retain age, gender, date of scan, date of birth, etc. which can be then used as regression covariates or to cluster patients together. Additional fields to remove can easily be added as an argument.

4. **Reliability** : a name with several typos or even a missing first name will still be anonymized correctly and linked to the correct demographics entry, by using fuzzy matching based on a normalized levenshstein distance for characters mistakes and a modified normalized words jaccard distance on words switching/missing.

5. **Robustness** : all other dicoms anonymizers look for the name only in the PatientName field. However, there are other hidden fields which can contain the patient's name, and this changes for each machine, each center, each parametrization of a scanner. This application will robustly find all fields containing the patient's name, and will try to remove it without touching the rest of the field (which can contain critical technical data), and it always checks that the name is cleared from the whole dicom before continuing. In case the patient name remains for whatever reason, the anonymizer will stop and show an error message so that you can take the appropriate steps.

To use this application, you must point to a folder where each subfolder is a patient dicom session. You can put multiple dicom sessions if it is for the same patient. This is because the script assumes one patient per folder (to look for name, else it would be too expensive to look at all dicoms files).

This app can work with a GUI or from commandline interchangeably, the same features and options will be available.

Another interesting feature is that you can build on top of an anonymized dataset, you can add new subjects. Indeed, each subject gets a unique id based on the name (by default) or by the order in the folder (by using the appropriate option), which makes it impossible to translate back to the original name from the id in any case.

TODO
---------
This application is provided as-is with the hopes it can be useful to someone else. It was not cleaned up beforehand more than what was sufficient for running as a standalone GUI application and packaged as a minimal size binary.

With that in mind, the following should be done to make this application a proper software:

* Make demographics file optional. It is already not really necessary in the code to anonymize dicoms, an empty csv with just the 'name' column should be enough. But it would be better if it was really optional and handled correctly.
* Refactor code: put functions above the code, delete the # [] useless comments coming from Jupyter export, move to functions some code snippets.
* Use requirements instead of integrating submodules in csg_fileutil_libs.
* Unit test with randomly generated dicoms.

LICENSE
-------------
CSG Dicoms Anonymizer was initially made by Stephen Larroque <LRQ3000> for the Coma Science Group - GIGA Consciousness - CHU de Liege, Belgium. The application is licensed under MIT License.


.. |LICENCE| image:: https://img.shields.io/pypi/l/csg_dicoms_anonymizer.svg
   :target: https://raw.githubusercontent.com/lrq3000/csg_dicoms_anonymizer/master/LICENCE
.. |PyPI-Status| image:: https://img.shields.io/pypi/v/csg_dicoms_anonymizer.svg
   :target: https://pypi.python.org/pypi/csg_dicoms_anonymizer
.. |PyPI-Versions| image:: https://img.shields.io/pypi/pyversions/csg_dicoms_anonymizer.svg
   :target: https://pypi.python.org/pypi/csg_dicoms_anonymizer
