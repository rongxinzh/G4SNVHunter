# NEWS for G4SNVHunter package

CHANGES IN VERSION 0.99.0
--------------
- Initial release of the G4SNVHunter package.

CHANGES IN VERSION 0.99.1
--------------
- Minor updates to documentation.

CHANGES IN VERSION 0.99.2
--------------
- Remove unnecessary project files from version control.

CHANGES IN VERSION 0.99.3
--------------
- Minor revisions were made based on Bioconductor reviewers' feedback.

CHANGES IN VERSION 0.99.4
--------------
- Minor updates

CHANGES IN VERSION 1.0.0
--------------
- Bioconductor auto-bumped release version without code changes.

CHANGES IN VERSION 1.1.0
--------------
- Bioconductor auto-bumped devel version to 1.1.0 without code changes.

CHANGES IN VERSION 1.1.1
------------------------

- Deprecated functions (still supported but scheduled for removal in a future release): 
  * checkSNV()
  * SNVImpactG4()         | replaced by G4VarImpact()
  * filterSNVImpact()     | replaced by filterVarImpact()
  * plotSNVImpact()       | replaced by plotVarImpact()
  * plotImpactSeq()       | replaced by plotImpactedG4()

- Updated functions:
  * G4HunterDetect()      | now includes global metadata in the returned object

- Removed functions:
  * validateG4HunterParams() | replaced by validateG4HunterDetectInputs()

- New functions:
  * exportG4()
  * exportMutG4()
  * filterVarImpact()
  * G4VarImpact()
  * loadVariant()
  * plotVarImpact()
  * plotImpactedG4()
  
CHANGES IN VERSION 1.1.2
------------------------
- Tweaked the vignette

CHANGES IN VERSION 1.1.3
------------------------
- Remove .o files

CHANGES IN VERSION 1.1.4
------------------------
- Attempting to resolve build issues introduced by the uniquely genius Taishan compilation environment
