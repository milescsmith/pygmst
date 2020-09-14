# 2020-09-14 - 0.4.19

## Fixes
* Found that file names with a `--` in them are treated by `probuild` as
  commandline arguments.  Pygmst now copies the sequence file to a temporary file
  (since symlinks don't always work) to avoid this problem.

## Additions
* Add CHANGELOG.md