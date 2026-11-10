### Short description
* > information may be structured by bullet points and may be used as commit message
* summary of core changes in science
* summary of core changes in technical aspects of the code
* do simulation results change (brief summary if so)


### Detailed description
#### Issue
* description of the issue / feature you want to address

#### Solution
* what has changed (science and technical aspects) in the code and where (process, module, routine/function)
  * please generally avoid adding chunks of the code to this MR description
  * instead provide a brief description of your changes along with `module-name:routine-name` (and optionally the line number)
* plus scientific justification, if applicable

#### Additional changes
* additional changes worth mentioning that do not directly belong to the issue / solution but are nevertheless connected to the issue in a way that they can be seen as part of the same feature

#### Follow-up or remaining issues
* add here a list of new (follow-up) or remaining issues, if applicable
* if suitable already create a gitlab issue provide the link to it

#### ICON-Land testing and docu of changes in simulation results
* [ ] I tested the latest feature-complete version of the branch (using the NAG & GPU compilers)
  * [ ] `exp.land_jsbach_R2B4_test`
  * [ ] `exp.land_jsbach_sitelevel_test`
  * [ ] `exp.land_quincy_canopy_R2B4_test`
  * [ ] `exp.land_iqj_sitelevel_test`
  * [ ] GPU and CPU simulatio results do NOT differ more than expected (as with previous model version)
* [ ] simulation results do NOT change
* [ ] simulation results do change

#### Associated merge requests
* ICON: add link (full url) to the release notes provided with the ICON MR
