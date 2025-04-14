## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Resubmission 2023-04-31
This is a resubmission. In this version I have:

* Removed the use of the package name BREADR in the description part of DESCRIPTION to avoid spell-checker errors. 

* Wrapped the example for saveSlices in \dontrun{} to avoid excess execution time. 

## Resubmission 2023-04-05

> If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

* Added a link to the preprint manuscript into the DESCRIPTION file. 

> \dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.
Does not seem necessary.
Please unwrap the examples if they are executable in < 5 sec, or replace
\dontrun{} with \donttest{}.

* Changed to donttest{}.

> You write information messages to the console that cannot be easily
suppressed. It is more R like to generate objects that can be used to
extract the information a user is interested in, and then print() that
object.
Instead of print()/cat() rather use message()/warning()  or
if(verbose)cat(..) (or maybe stop()) if you really have to write text to
the console.
(except for print, summary, interactive functions)

* Changed the output from test_degree to print to screen only if verbose is set to TRUE. 

## Resubmission 2023-04-06

> Please write references in the description of the DESCRIPTION file in
the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: authors (year) <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

* Added author and year - used et al. as four authors. 

> You still write information messages to the console that cannot be
easily suppressed. It is more R like to generate objects that can be
used to extract the information a user is interested in, and then
print() that object.
Instead of print()/cat() rather use message()/warning()  or
if(verbose)cat(..) (or maybe stop()) if you really have to write text to
the console.
(except for print, summary, interactive functions)
e.g.: R/plotLOAF.R

* Changed all output to be wrapped in verbose and added verbose as parameter. 

## Resubmission 2023-04-07

> Please also add the verbose argument to all functions with cat() as you correctly did for R/plotLOAF.R

> -> e.g.: R/processEigenstrat.R, R/saveSLICES.R...

* All cat functions wrapped in verbose if condition. 

## Patch 2023-04-12

* Found an error in read functions that causing error in PMR. Fixed. 

## Patch 2024-09-09

* Fixed processEigenstrat so that works on windows machines. 

## Resubmission 2024-09-09

Mistakenly left repdev checking in folder - now removed. 

## Patch 2025-04-14

Submitting paper to JOSS. Reviewer could only install if R4.4, so changed to this from their recommendation. 