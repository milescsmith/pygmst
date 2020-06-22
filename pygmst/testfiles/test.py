# import sys
# import click
# from logger import Logipy

# verbose = False


# @click.command()
# @click.option("--verbose", is_flag=True, default=False, show_default=True)
# @click.option("--input", type=str, required=True)
# @click.help_option(show_default=True)
# @Log(verbose)
# def main(input, verbose):
#     print(input)


# if __name__ == "__main__":
#     if "--verbose" in sys.argv:
#         verbose = True
#     sys.exit(main())  # pragma: no cover
