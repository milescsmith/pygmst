# [0.5.0] - 2021-05-11

## Changed
- Switched from using `click` to `typer` for parsing the cli
- Replace usage of os module with pathlib
- Now using poetry for installation and dependency handling so that we are
  PEP517 and PEP633 compliant
- Eliminated unnecessary `main()` function and now passing directly to `gmst()`
  (since that is *all* `main()` was doing)

## Fixed
- gmst() now checks that the file size for the sequence file isn't 0

# [0.4.19] - 2020-09-14

## Fixes
- Found that file names with a `--` in them are treated by `probuild` as
  commandline arguments.  Pygmst now copies the sequence file to a temporary file
  (since symlinks don't always work) to avoid this problem.

## Additions
* Add CHANGELOG.md