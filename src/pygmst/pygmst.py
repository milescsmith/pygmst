#!/usr/bin/env python

import logging
import re
from tempfile import NamedTemporaryFile
from shutil import copyfile
from subprocess import check_output, CalledProcessError
from typing import Dict, List, Optional, Tuple
from pathlib import Path


import pyfaidx
from Bio.SeqUtils import GC as getGC
from pkg_resources import resource_filename
from sortedcontainers import SortedDict
import typer
from enum import Enum, IntEnum

from pygmst import __version__


app = typer.Typer(name="pygmst", help="Python translation of GMST")


class StrandOption(str, Enum):
    direct = "direct"
    reverse = "reverse"
    both = "both"


class TranslationTableOption(str, Enum):
    eleven = "11"
    four = "4"
    one = "1"


class MotifOption(str, Enum):
    zero = "0"
    one = "1"


class GibbsOption(str, Enum):
    one = "1"
    three = "3"


seq = "sequence"
start_prefix = "startseq."
gibbs_prefix = "gibbs_out."
mod_prefix = "itr_"
mod_suffix = ".mod"
hmmout_prefix = "itr_"
hmmout_suffix = ".lst"
out_name = "GeneMark_hmm.mod"
fnn_out = ""
faa_out = ""

meta_out = "initial.meta.lst"
gc_out = f"{meta_out}.feature"
verbose = False


def version_callback(value: bool):
    """Prints the version of the package."""
    if value:
        print(f"pygmst version: {__version__}")
        raise typer.Exit()


def setup_logging(name: Optional[str] = None):
    if name:
        logger = logging.getLogger(name)
    else:
        logger = logging.getLogger(__name__)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    if name:
        fh = logging.FileHandler(filename=name)
    else:
        fh = logging.FileHandler(filename=f"{__name__}.log")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    st = logging.StreamHandler()
    st.setLevel(logging.INFO)
    st.setFormatter(formatter)
    logger.addHandler(st)


@app.command(name="")
def gmst(
    seqfile: str = typer.Argument(...),
    output: Optional[str] = typer.Option(
        None,
        "--output",
        "-o",
        help=(
            "output file with predicted gene coordinates by GeneMark.hmm and"
            "species parameters derived by GeneMarkS-T. If not provided, the"
            "file name will be taken from the input file. GeneMark.hmm can be"
            "executed independently after finishing GeneMarkS training."
            "This method may be the preferable option in some situations, as"
            "it provides accesses to GeneMark.hmm options."
        ),
    ),
    outputformat: Optional[str] = typer.Option(
        "LST",
        "--format",
        "--outputformat",
        help="output coordinates of predicted genes in LST or GTF format.",
    ),
    fnn: bool = typer.Option(
        False, help="create file with nucleotide sequence of predicted genes"
    ),
    faa: bool = typer.Option(
        False, help="create file with protein sequence of predicted genes"
    ),
    clean: bool = typer.Option(True, help="delete all temporary files"),
    bins: int = typer.Option(
        0,
        min=1,
        max=4,
        clamp=True,
        help="number of clusters for inhomogeneous genome. Use 0 for automatic clustering",
    ),
    prok: bool = typer.Option(
        False,
        help="to run program on prokaryotic transcripts (this option is the same as: --bins 1 --filter 0 --order 2 --order_non 2 --gcode 11 -width 6 --prestart 40 --fixmotif 0)",
    ),
    filterseq: int = typer.Option(
        1, "--filter", help="keep at most one prediction per sequence"
    ),
    strand: StrandOption = typer.Option(
        StrandOption.both, help="sequence strand to predict genes in"
    ),
    order: int = typer.Option(4, min=1, max=100, clamp=True, help="markov chain order"),
    order_non: int = typer.Option(2, help="order for non-coding parameters"),
    gcode: TranslationTableOption = typer.Option(
        TranslationTableOption.one,
        help="genetic code.  See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for codes.  Currently, only tables 1, 4, and 11 are allowed.",
    ),
    motifopt: MotifOption = typer.Option(
        MotifOption.one,
        "--motif",
        help="iterative search for a sequence motif associated with CDS start",
    ),
    width: int = typer.Option(12, min=3, max=100, help="motif width"),
    prestart: int = typer.Option(
        6,
        min=0,
        max=100,
        help="length of sequence upstream of translation initiation site that presumably includes the motif",
    ),
    fixmotif: bool = typer.Option(
        True,
        help="the motif is located at a fixed position with regard to the start motif could overlap start codon. if this option is on, it changes the meaning of the --prestart option which in this case will define the distance from start codon to motif start",
    ),
    offover: bool = typer.Option(True, help="prohibits gene overlap"),
    par: Optional[str] = typer.Option(
        None,
        help="custom parameters for GeneMarkS (default is selected based on gcode value: 'par_<gcode>.default' )",
    ),
    gibbs: GibbsOption = typer.Option(
        GibbsOption.three,
        help="version of Gibbs sampler software (default: 3 supported versions: 1 and 3 )",
    ),
    test: bool = typer.Option(False, help="installation test"),
    identity: float = typer.Option(
        0.99,
        min=0.0,
        max=1.0,
        help="identify level assigned for termination of iterations",
    ),
    maxitr: int = typer.Option(
        10,
        min=1,
        help="maximum number of iterations (default: 10 supported in range: >= 1)",
    ),
    verbose: int = typer.Option(0, "--verbose", "-v", count=True),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    """Front end to GeneMark S-T.  Predict open reading frames from raw sequences

    \b
    Parameters
    ----------
    seqfile:
        FASTA containing sequences to use for prediction

    """
    setup_logging("pygmst")

    if verbose == 1:
        logging.basicConfig(filename="pygmst.log", filemode="w", level=logging.WARN)
    elif verbose == 2:
        logging.basicConfig(filename="pygmst.log", filemode="w", level=logging.INFO)
    elif verbose == 3:
        logging.basicConfig(filename="pygmst.log", filemode="w", level=logging.DEBUG)

    motifopt = int(motifopt.value)
    gcode = int(gcode.value)
    gibbs = int(gibbs.value)

    motif = bool(motifopt)  # for compatibility
    fixmotif = True if fixmotif == 1 else False  # for compatibility

    seqfile = Path(seqfile)
    if output is None:
        if outputformat == "LST":
            output = seqfile.with_suffix(".lst")
        elif outputformat == "GFF":
            output = seqfile.with_suffix(".gff")
        if fnn:
            fnn_out = seqfile.with_suffix(".fnn")
        if faa:
            faa_out = seqfile.with_suffix(".faa")
    else:
        if fnn:
            fnn_out = seqfile.with_suffix(".fnn")
        if faa:
            faa_out = seqfile.with_suffix(".faa")
    if prok:
        bins = 1
        filterseq = 0
        order = 2
        order_non = 2
        offover = False
        gcode = 11
        fixmotif = False
        prestart = 40
        width = 6
    if par is None:
        par = resource_filename("pygmst", f"genemark/par_{gcode}.default")

    # with tempfile.TemporaryDirectory() as tmpdir:
    # tmpdir = tempfile.mkdtemp()
    # out_name = f"{tmpdir}/{out_name}"
    # GeneMark.hmm gene finding program <gmhmmp>; version 2.14
    hmm = resource_filename("pygmst", "genemark/gmhmmp")

    # sequence parsing tool <probuild>; version 2.11c
    build = resource_filename("pygmst", "genemark/probuild")

    # gibbs sampler - from http://bayesweb.wadsworth.org/gibbs/gibbs.html ; version 3.10.001  Aug 12 2009
    # this doesn't appear to be available from that source?
    # possible replacement at https://github.com/Etschbeijer/GibbsSampler ?
    gibbs3 = resource_filename("pygmst", "genemark/Gibbs3")

    # use <probuild> with GeneMarkS parameter file <$par>
    build = f"{build} --par"

    # set options for <gmhmmp>

    # switch gene overlap off in GeneMark.hmm; for eukaryotic intron-less genomes
    if offover:
        hmm = f"{hmm} -p 0"

    # set strand to predict
    if strand == "direct":
        hmm = f"{hmm} -s d"
    elif strand == "reverse":
        hmm = f"{hmm} -s r"

    ## tmp solution: get sequence size, get minimum sequence size from --par <file>
    ## compare, skip iterations if short

    # this should end up as something like
    # 'probuild --par par_1.default --clean_join sequence --seq test.fa

    # okay, so probuild cannot tolerate a file with two dashes in the name
    # I believe it thinks you are trying to pass an argument since it gives the
    # generic "error in command line".  So, create a symlink with a non-offending
    # name and use that

    if Path(seqfile).stat().st_size == 0:
        raise ValueError(f"The input sequence {Path(seqfile).resolve()}is empty")

    sinnombre = copyfile(src=seqfile, dst=NamedTemporaryFile().name)
    probuild_cmd = f"{build} {par} --clean_join {seq} --seq {sinnombre}"

    check_output(probuild_cmd.split())

    with open(seqfile, "r") as _:
        sequence_size = len(_.read())

    command = f"grep MIN_SEQ_SIZE {par}"
    result = check_output(args=command.split())
    minimum_sequence_size = int(
        re.findall(
            pattern=r"\s*--MIN_SEQ_SIZE\s+(\d+)",
            string=str(result, "utf=8"),
        )[0]
    )

    do_iterations = 1

    if sequence_size < minimum_sequence_size:
        do_iterations = 0

    if do_iterations > 0:
        # ------------------------------------------------
        # clustering
        # initial prediction using MetaGeneMark model

        # get GC of whole sequence for each sequence"
        gc_out = f"{meta_out}.feature"

        # form of! 'probuild --par par_1.default --stat_fasta --seq test.fa > initial.meta.list.feature'
        gc_cmd = f"{build} {par} --stat_fasta --seq {sinnombre}"
        with open(gc_out, "w+") as gc_capture:
            logging.debug(gc_cmd)
            gc_capture.write(str(check_output(gc_cmd.split()), "utf-8"))

        # determine bin number and range
        bin_num, cutoffs, seq_GC = cluster(
            feature_f=gc_out, clusters=bins, min_length=10000
        )

        logging.info(f"bin number = {bin_num}\n")
        logging.info(f"GC range = {''.join([str(_) for _ in cutoffs])}\n")

        # training
        models: List[str] = list()
        seqs: List[str] = list()
        # my %handles; # file handles for n bins.
        if bin_num == 1:
            logging.debug(
                f"train(input_seq={sinnombre}, "
                f"fseq={seq}, "
                f"motif={motif},"
                f"fixmotif={fixmotif}, "
                f"order={order}, "
                f"order_non={order_non}, "
                f"start_prefix={start_prefix}, "
                f"gibbs_prefix={gibbs_prefix}, "
                f"prestart={prestart}, "
                f"width={width}, "
                f"build_cmd={build}, "
                f"hmm_cmd={hmm}, "
                f"par={par}, "
                f"maxitr={maxitr}, "
                f"identity={identity}, "
                f"gibbs3={gibbs3})"
            )
            final_model = train(
                input_seq=seqfile,
                seq=seq,
                motif=motif,
                fixmotif=fixmotif,
                order=order,
                order_non=order_non,
                start_prefix=start_prefix,
                gibbs_prefix=gibbs_prefix,
                prestart=prestart,
                width=width,
                build_cmd=build,
                hmm_cmd=hmm,
                par=par,
                maxitr=maxitr,
                identity=identity,
                gibbs3=gibbs3,
                # tmpdir=tmpdir,
            )
            logging.info(f"final_model: {final_model}")

            # make a GC file for the input file
            # read input sequences
            try:
                FA = pyfaidx.Fasta(seqfile, duplicate_action="longest")
                seq_GC_entries = [_.lstrip(">") for _ in seq_GC.keys()]
                for read in FA:
                    if read.long_name not in seq_GC_entries:
                        seq_GC[f">{read.long_name}"] = int(getGC(read[:].seq))
            except IOError as e:
                logging.critical(f"{e}\nCannot open {seqfile}")

        else:
            # create sequence file for each bin
            logging.info("Binning the input sequence")

            # for i in range(0, bin_num):
            #     with open(file=f"{tmpdir}/seq_bin_{i}", mode="w") as fh:
            #         seqs.extend([f"{tmpdir}/seq_bin_{i}"])
            #         handles[i] = f"{tmpdir}/seq_bin_{i}"
            #         logging.info(f"Created {tmpdir}/seq_bin_{i}")

            # read input sequences
            FA = pyfaidx.Fasta(seqfile, duplicate_action="longest")
            seq_GC_entries = [_.lstrip(">") for _ in seq_GC.keys()]
            seq_bins: Dict[int, List[str]] = {
                k: [] for k in (_ for _ in range(0, bin_num))
            }  # NEW! we are going to store the names of each read in a dictionary, which has a list member for each bin
            for read in FA:  # for each record
                if (
                    read.long_name not in seq_GC_entries
                ):  # if this has not already had the %GC determined
                    seq_GC[read.long_name] = int(getGC(read[:].seq))  # Do it now

                    # decide which bin the sequence belongs to, add the name to the bin dictionary
                if (
                    bin_num == 2
                ):  # if there are two bins, we can really only be above or below a cutoff
                    if seq_GC[read.long_name] <= cutoffs[0]:
                        seq_bins[0].extend([read.long_name])
                    else:
                        seq_bins[1].extend([read.long_name])
                else:  # else, does the read have low, medium, or high %GC?
                    if seq_GC[read.long_name] <= cutoffs[0]:
                        seq_bins[0].extend([read.long_name])
                    elif seq_GC[read.long_name] <= cutoffs[1]:
                        seq_bins[1].extend([read.long_name])
                    else:
                        seq_bins[2].extend([read.long_name])
                # logging.info(f"Placed {read_header} into bin {bin}")

            # write all sequences to a bin file at once
            try:
                for i in range(0, bin_num):
                    with open(file=f"seq_bin_{i}", mode="w") as fh:
                        for j in seq_bins[i]:
                            fh.writelines(
                                f"{FA[j].long_name}\t[gc={seq_GC[j]}]\n{FA[j][:].seq}\n"  # modify the read.long_name to include the %GC
                            )
                    logging.info(f"Created seq_bin_{i}")

                seqs = [f"seq_bin_{i}" for i in range(0, bin_num)]
                # handles = SortedDict({k:f"{tmpdir}/seq_bin_{k}" for k in range(0, bin_num)})
            except IOError as e:
                logging.critical(f"{e}\nCannot open {seqfile}")
                print(f"Cannot open {seqfile}")

            # train
            for i in range(bin_num):
                logging.debug(f"training on {seqs[i]}")
                current_model = train(
                    input_seq=seqs[i],
                    seq=seq,
                    motif=motif,
                    fixmotif=fixmotif,
                    order=order,
                    order_non=order_non,
                    start_prefix=start_prefix,
                    gibbs_prefix=gibbs_prefix,
                    prestart=prestart,
                    width=width,
                    build_cmd=build,  # probuild --par par_1.default
                    hmm_cmd=hmm,
                    par=par,
                    maxitr=maxitr,
                    identity=identity,
                    gibbs3=gibbs3,
                    bin_num=i,
                    # tmpdir=tmpdir,
                )
                logging.info(f"train() returned: {current_model}")
                if current_model:
                    models.extend([current_model])

            logging.debug(
                f"combine {len(models)} individual models to make the final model file"
            )
            model_names = "\n".join(models)
            logging.debug(f"those models are: {model_names}")
            final_model = combineModels(
                mod=models,
                cut_offs=cutoffs,
            )  # tmpdir=tmpdir)
            logging.info(f"combineModels() returned: {final_model}")

        check_output(f"cp {final_model} {out_name}".split())
        logging.info(f"Ran 'cp {final_model} {out_name}' but not sure why?")

    if outputformat == "GFF":
        format_gmhmmp = " -f G "
    else:
        format_gmhmmp = ""

    if filterseq == 1:
        hmm += " -b "
    if faa:
        hmm += f" -A {faa_out} "
    if fnn:
        hmm += f" -D {fnn_out} "
    # $hmm .= " -a $faa_out " if $faa;
    # $hmm .= " -d $fnn_out " if $fnn;

    if do_iterations:
        if motif:
            command = f"{hmm} -r -m {out_name} -o {output} {format_gmhmmp} {seqfile}"
            result = check_output(command.split())
            logging.debug(command)
            logging.debug(f"{str(result, 'utf-8')}")
        else:
            # no moitf option specified
            command = f"{hmm} -m {out_name} -o {output} {format_gmhmmp} {seqfile}"
            result = check_output(command.split())
            logging.debug(command)
            logging.debug(f"{str(result, 'utf-8')}")
    else:
        meta_model = resource_filename("pygmst", "genemark/MetaGeneMark_v1.mod")
        # no iterations - use heuristic only
        command = f"{hmm} -m {meta_model} -o {output} {format_gmhmmp} {seqfile}"
        result = check_output(command.split())
        logging.debug(command)
        logging.debug(f"{str(result, 'utf-8')}")

    print(f"wrote final results to {output}")


def train(
    input_seq: str,  # test.fa
    seq: str,  # sequence
    motif: bool,  # True
    fixmotif: bool,  # True
    order: int,  # 4
    order_non: int,  # 2
    start_prefix: str,  # "startseq."
    gibbs_prefix: str,  # "itr_"
    prestart: int,  # 6
    width: int,  # 12
    build_cmd: str,  # "probuild --par"
    hmm_cmd: str,  # "gmhmmp -s d"
    par: str,  # "par_1.default"
    maxitr: int,  # 10
    identity: float,  # 0.99
    gibbs3: str,
    bin_num: int = 0,
    # tmpdir: Optional[str] = None,
) -> str:
    logging.info("Beginning training")
    # ------------------------------------------------
    # prepare sequence
    # build_cmd should have the form of! `probuild --par par_1.default`
    # so this makes the below
    # `probuild --par par_1.default --clean_join sequence --seq test.fa`
    # which seems kind of redundant, but what do I know

    if Path(input_seq).stat().st_size == 0:
        raise ValueError("the input sequence is empty")

    command = f"{build_cmd} {par} --clean_join {seq} --seq {input_seq}"
    result = check_output(command.split())
    logging.debug(command)
    logging.debug(f"{str(result, 'utf-8')}")

    # tmp solution: get sequence size, get minimum sequence size from --par <file>
    # compare, skip iterations if short

    with open(seq, "r") as _:
        sequence_size = len(_.read())

    min_size_find_cmd = check_output(args=f"grep MIN_SEQ_SIZE {par}".split())
    minimum_sequence_size = int(
        re.findall(
            pattern=r"\s*--MIN_SEQ_SIZE\s+(\d+)",
            string=str(min_size_find_cmd, "utf=8"),
        )[0]
    )

    do_iterations = True
    itr = 0

    logging.info(f"minimum_sequence_size: {minimum_sequence_size}")
    logging.info(f"sequence_size: {sequence_size}")
    if sequence_size < minimum_sequence_size:
        # form of! `probuild --par par_1.default --clean_join sequence --seq test.fa --MIN_CONTIG_SIZE 0 --GAP_FILLER
        command = f"{build_cmd} {par} --clean_join {seq} --seq {input_seq} --MIN_CONTIG_SIZE 0 --GAP_FILLER"
        result = check_output(command.split())
        logging.debug(command)
        logging.debug(f"{str(result, 'utf-8')}")

        init_mod = resource_filename("pygmst", "genemark/MetaGeneMark_v1.mod")
        next_item = f"bin_{bin_num}{hmmout_suffix}"
        command = f"{hmm_cmd} -m {init_mod} -o {next_item} {seq}"
        result = check_output(command.split())
        logging.debug(command)
        logging.debug(f"{str(result, 'utf-8')}")

        mod = f"bin_{bin_num}{mod_suffix}"
        # $command = "$build            --mkmod $mod  --seq $seq  --geneset $next       --ORDM $order  --order_non $order_non  --revcomp_non 1";
        command = f"{build_cmd} {par} --mkmod {mod} --seq {seq} --geneset {next_item} --ORDM {order} --order_non {order_non} --revcomp_non 1"

        if motif and not fixmotif:
            command = (
                f"{command} --pre_start {start_prefix} --PRE_START_WIDTH {prestart}"
            )
        elif motif and fixmotif:
            command = (
                f"{command} --fixmotif --PRE_START_WIDTH {prestart} --width {width}"
            )

        try:
            results = check_output(command.split())
            logging.debug(command)
            logging.debug(f"{str(result, 'utf-8')}")
        except CalledProcessError as e:
            logging.debug(e)
        return mod  # this may end up returning an empty file because the bin{bin_num}.lst file is empty or has insufficient number of entries.
        # this, in turn, fouls up the

    # Log("do_iterations = $do_iterations\n")

    # ------------------------------------------------
    # run initial prediction
    next_item = f"{hmmout_prefix}{itr}_bin_{bin_num}{hmmout_suffix}"
    print("run initial prediction")

    # form of! `gmhmmp -s d sequence -m MetaGeneMark_v1.mod -o itr_{itr}.lst`
    mod = resource_filename("pygmst", "genemark/MetaGeneMark_v1.mod")
    command = f"{hmm_cmd} -m {mod} -o {next_item} {seq}"
    result = check_output(command.split())
    logging.debug(command)
    logging.debug(f"{str(result, 'utf-8')}")

    # enter iterations loop
    if do_iterations:
        logging.info("entering training iteration loop")
    while do_iterations:
        itr += 1
        mod = f"{mod_prefix}{itr}_bin_{bin_num}{mod_suffix}"
        logging.info(f"iteration: {itr}, mod: {mod}")

        if motif and not fixmotif:
            start_seq = f"{start_prefix}{itr}"
            gibbs_out = f"{gibbs_prefix}{itr}"

        command = f"{build_cmd} {par} --mkmod {mod} --seq {seq} --geneset {next_item} --ORDM {order} --order_non {order_non} --revcomp_non 1"

        if motif and not fixmotif:
            command = f"{command} --pre_start {start_seq} --PRE_START_WIDTH {prestart}"
        elif motif and fixmotif:
            command = (
                f"{command} --fixmotif --PRE_START_WIDTH {prestart} --width {width}"
            )
            # command += f" --fixmotif --PRE_START_WIDTH {prestart} --width {width} --log {logfile}"

        logging.info(f"build model: {mod} for iteration: {itr}")
        results = check_output(command.split())
        logging.debug(command)
        logging.debug(f"{str(result, 'utf-8')}")

        if (
            motif and not fixmotif
        ):  # given that we only have Gibbs3, don't *really* have an option
            # if $gibbs_version == 1:
            #     &RunSystem( "$gibbs $start_seq $width -n > $gibbs_out", "run gibbs sampler\n" )
            # elif $gibbs_version == 3:
            #     &RunSystem( "$gibbs3 $start_seq $width -o $gibbs_out -F -Z  -n -r -y -x -m -s 1 -w 0.01", "run gibbs3 sampler\n" )
            logging.info("run gibbs3 sampler")
            # form of! 'Gibbs3 startseq.{itr} 12 -o gibbs_out.{itr} -F -Z -n -r -y -x -m -s 1 -w 0.01"
            command = f"{gibbs3} {start_seq} {width} -o {gibbs_out} -F -Z -n -r -y -x -m -s 1 -w 0.01"
            result = check_output(command.split())
            logging.debug(command)
            logging.debug(f"{str(result, 'utf-8')}")

            logging.info("make prestart model")
            # run([build_cmd, "--gibbs", gibbs_out, "--mod", mod, "--seq", start_seq, "--log", logfile])
            # form of! `probuild --par par_1.default --gibs gibbs_out.{itr} --mod itr_{itr}.mod --seq startseq.{itr}
            command = (
                f"{build_cmd} {par} --gibbs {gibbs_out} --mod {mod} --seq {start_seq}"
            )
            result = check_output(command.split())
            logging.debug(command)
            logging.debug(f"{str(result, 'utf-8')}")

        prev = next_item
        # next_item = itr_{itr}.lst
        next_item = f"{hmmout_prefix}{itr}_bin_{bin_num}{hmmout_suffix}"

        # form of! `gmhmmp -s d sequence -m itr_{itr}.mod -o itr_{itr}.lst
        command = f"{hmm_cmd} {seq} -m {mod} -o {next_item}"

        # form of! `gmhmmp -s d sequence -m itr_{itr}.mod -o itr_{itr}.lst -r
        if motif:
            command += " -r"

        logging.info(f"prediction, iteration: {itr}")
        result = check_output(command.split())
        logging.debug(command)
        logging.debug(f"{str(result, 'utf-8')}")

        # `probuild --par par_1.default --compare --source itr_{itr}.lst --target itr_{itr-1}.lst
        command = f"{build_cmd} {par} --compare --source {next_item} --target {prev}"
        # &Log( "compare:\n" . $command . "\n" );

        results = check_output(command.split())
        logging.debug(command)
        logging.debug(f"{str(result, 'utf-8')}")

        diff = str(results, "utf-8").strip("\n")
        logging.info(f"iteration {itr}, bin {bin_num} difference: {diff}")
        # &Log( "compare $prev and $next_item: $diff\n" );

        if float(diff) >= identity:
            logging.info(f"Stopped iterations on identity: {diff}")
            do_iterations = False
        if itr == maxitr:
            logging.info(f"Stopped iterations on maximum number: {maxitr}")
            return mod
    logging.info(f"return value of train: {mod}")
    return mod


def cluster(
    feature_f: str, clusters: int, min_length: int = 10000
) -> Tuple[int, List[int], Dict[str, int]]:  # $gc_out, $bins
    logging.debug("Beginning clustering")
    gc_hash: Dict[int, int] = dict()
    cut_off_points: List[int] = list()
    num_of_seq = 0
    total_length = 0
    header_to_cod_GC = dict()

    # with feature_f as GC:
    # read in probuild output, line by line.  Should be fasta input.
    with open(feature_f, "r") as GC:
        for line in GC:
            # if the line is a fasta header in the form of '>(Reference sequence name)\t(number) (number)
            # if (text := re.search(pattern="^>(.*?)\t(\d+)\s+(\d+)", string=line)): # switch to this?  only support python>=3.8?
            text = re.search(pattern=r"^>(.*?)\t(\d+)\s+(\d+)", string=line)
            if text:
                # logging.debug(f"text.group(0) = {text.group(0)}")
                header = text.group(1)  # Reference name
                # logging.debug(f"header = {header}")
                length = int(text.group(2))  # length of sequence?
                # logging.debug(f"header = {length}")
                seqGC = int(text.group(3))  # must be GC percentage
                # logging.debug(f"header = {seqGC}")

                header_to_cod_GC[header] = seqGC
                num_of_seq += 1
                total_length += length
                if seqGC in gc_hash:
                    gc_hash[seqGC] += length
                else:
                    gc_hash[seqGC] = length

    sorted_GC = SortedDict(gc_hash)  # sort the gc_hash dictionary by keys
    min_GC = sorted_GC.keys()[0]
    max_GC = sorted_GC.keys()[-1]
    print(f"min_GC={min_GC} max_GC={max_GC} total_seq_length={total_length}\n")

    previous = 0
    for key in sorted_GC:
        gc_hash[key] += previous

        if previous < total_length / 3 and gc_hash[key] >= total_length / 3:
            one_third = key
            logging.debug(f"({one_third})->({gc_hash[one_third]})\n")
        if previous < total_length / 3 * 2 and gc_hash[key] >= total_length / 3 * 2:
            two_third = key
            logging.debug(f"({two_third})->({gc_hash[two_third]})\n")
        if previous < total_length / 2 and gc_hash[key] >= total_length / 2:
            one_half = key
            logging.debug(f"({one_half})->({gc_hash[one_half]})\n")
        previous = gc_hash[key]

    if clusters == 0:
        # cluster number is not specified by user
        # automatically choose cluster number.
        if two_third - one_third > 3:
            clusters = 3
        else:
            clusters = 1
    if clusters == 3:
        if (
            (two_third - one_third) < 1
            or (max_GC - two_third) < 1
            or (one_third - min_GC) < 1
        ):
            # &Log( "Total number of sequences is not enough for training in 3 clusters!\n" )
            clusters = 1
        else:
            if gc_hash[one_third] > min_length:
                cut_off_points.extend([min_GC, one_third, two_third, max_GC])
            else:
                # &Log( "Total length of sequences is not enough for training in 3 clusters!\n" )
                clusters = 2

    if clusters == 2:
        if gc_hash[one_half] > min_length:
            cut_off_points.extend([min_GC, one_half, max_GC])
        else:
            # &Log( "Total length of sequences is not enough for training in 2 clusters!\n" )
            pass

    if clusters == 1:
        cut_off_points.extend([min_GC, max_GC])
    return clusters, cut_off_points, header_to_cod_GC


def combineModels(
    mod: List[str],
    cut_offs: List[int],
    minGC: int = 30,
    maxGC: int = 70,
    # tmpdir: str = ".",
) -> str:
    """
    concatenate model files into one in MGM format:
    starts with "__GC" and ends with "end"
    """
    logging.debug("combining models")
    model_names = "\n".join(mod)
    logging.debug(f"models for combining: {model_names}")

    # change the min and max GC value of cut_offs to the minGC and maxGC used by gmhmmp
    cut_offs[0] = minGC
    cut_offs[-1] = maxGC

    b = 1
    data = str()
    final_model = "final_model.mod"
    with open(file=final_model, mode="w+") as model:
        for i in range(minGC, maxGC + 1):
            model.write(f"__GC{i}\t\n")
            if i == cut_offs[b] or i == maxGC:
                logging.debug(f"combineModels processing {mod[b-1]}")
                if Path(mod[b - 1]).exists():
                    with open(file=mod[b - 1], mode="r") as fh:
                        for line in fh:
                            if "NAME" in line:
                                line = line.strip("\n")
                                line = f"{line}_GC<={i}\n"
                            data += line
                    model.write(data)
                    model.write("end \t\n\n")
                else:
                    logging.debug(f"combineModels could not find {mod[b-1]}")
                if b > len(mod):
                    break
                b += 1
    return final_model


if __name__ == "__main__":
    app()
