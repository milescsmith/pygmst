__version__ = "0.1.0"

import sys
import os
import click
import re
import logging
from click_option_group import optgroup
from typing import Optional, List, Tuple, Dict
from subprocess import run, PIPE, STDOUT
import pyfaidx
from Bio.SeqUtils import GC as getGC
from sortedcontainers import SortedDict
import tempfile
from pkgutil import get_data
from pkg_resources import resource_filename

# MetaGeneMark model for initial prediction
# ------------------------------------------------

# TODO: redefine these as needed within the code
# global seq
# global start_prefix
# global gibbs_prefix
# global mod_prefix
# global mod_suffix
# global hmmout_prefix
# global hmmout_suffix
# global out_name
# global fnn_out
# global faa_out
# global working_path

# global meta_out
# global gc_out

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


@click.command()
@optgroup.group("Output options", help="output is in current working directory")
@optgroup.option(
    "--output",
    "-o",
    type=str,
    help="output file with predicted gene coordinates by GeneMarh.hmm and species parameters derived by GeneMarkS-T. If not provided, the file name will be taken from the input file. GeneMark.hmm can be executed independently after finishing GeneMarkS training.This method may be the preferable option in some situations, as it provides accesses to GeneMarh.hmm options.",
    required=False,
)
@optgroup.option(
    "--format",
    "outputformat",
    type=str,
    help="output coordinates of predicted genes in LST or GTF format.",
    default="LST",
    show_default=True,
)
@optgroup.option(
    "--fnn",
    is_flag=True,
    help="create file with nucleotide sequence of predicted genes",
    default=False,
    show_default=True,
)
@optgroup.option(
    "--faa",
    is_flag=True,
    help="create file with protein sequence of predicted genes",
    default=False,
    show_default=True,
)
@optgroup.option(
    "--clean",
    is_flag=True,
    help="delete all temporary files",
    default=True,
    show_default=True,
)
# Run options:
@optgroup.group("Run options")
@optgroup.option(
    "--bins",
    type=click.IntRange(1, 4, clamp=True),
    help="number of clusters for inhomogeneous genome. Use 0 for automatic clustering",
    default=0,
    show_default=True,
)
@optgroup.option(
    "--filter",
    "filterseq",
    type=int,
    help="keep at most one prediction per sequence",
    show_default=True,
    default=True,
)
@optgroup.option(
    "--strand",
    type=click.Choice(["direct", "reverse", "both"]),
    help="sequence strand to predict genes in",
    default="both",
    show_default=True,
)
@optgroup.option(
    "--order",
    type=click.IntRange(1, 100, clamp=True),
    help="markov chain order",
    default=4,
    show_default=True,
)
@optgroup.option(
    "--order_non",
    type=int,
    help="order for non-coding parameters",
    default=2,
    show_default=True,
)
@optgroup.option(
    "--gcode",
    type=click.Choice(["11", "4", "1"]),
    help="genetic code.  See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for codes.  Currently, only tables 1, 4, and 11 are allowed.",
    default="1",
    show_default=True,
)
@optgroup.option(
    "--motif",
    "motifopt",
    type=click.Choice(["0", "1"]),
    help="iterative search for a sequence motif associated with CDS start",
    show_default=True,
    default="1",
)
@optgroup.option(
    "--width",
    type=click.IntRange(3, 100),
    help="motif width",
    show_default=True,
    default=12,
)
@optgroup.option(
    "--prestart",
    type=click.IntRange(0, 100),
    help="length of sequence upstream of translation initiation site that presumably includes the motif",
    show_default=True,
    default=6,
)
@optgroup.option(
    "--fixmotif",
    type=bool,
    is_flag=True,
    help="the motif is located at a fixed position with regard to the start motif could overlap start codon. if this option is on, it changes the meaning of the --prestart option which in this case will define the distance from start codon to motif start",
    show_default=True,
    default=True,
)
@optgroup.option(
    "--offover",
    type=bool,
    is_flag=True,
    default=True,
    help="prohibits gene overlap",
    show_default=True,
)
@optgroup.group("Combined output and run options")
@optgroup.option(
    "--prok",
    is_flag=True,
    default=False,
    help="to run program on prokaryotic transcripts (this option is the same as: --bins 1 --filter 0 --order 2 --order_non 2 --gcode 11 -width 6 --prestart 40 --fixmotif 0)",
    show_default=True,
)
@optgroup.group("Test/developer options")
@optgroup.option(
    "--par",
    type=str,
    help="custom parameters for GeneMarkS (default is selected based on gcode value: 'par_<gcode>.default' )",
    show_default=True,
)
@optgroup.option(
    "--gibbs",
    type=click.Choice(["1", "3"]),
    default="3",
    help="version of Gibbs sampler software (default: 3 supported versions: 1 and 3 )",
    show_default=True,
)
@optgroup.option("--test", is_flag=True, default=False, help="installation test")
@optgroup.option(
    "--identity",
    type=click.FloatRange(min=0, max=1, clamp=False),
    default=0.99,
    help="identity level assigned for termination of iterations",
    show_default=True,
)
@optgroup.option(
    "--maxitr",
    type=int,
    help="maximum number of iterations (default: 10 supported in range: >= 1)",
    show_default=True,
    default=10,
)
@optgroup.option("-v", "--verbose", count=True)
@optgroup.option("--version", is_flag=True, default=False, show_default=True)
@click.argument("seqfile", type=str)
@click.help_option(show_default=True)
# @Log(verbose)
def main(
    seqfile: str,
    output: str,
    outputformat: Optional[str] = None,
    fnn: bool = False,
    faa: bool = False,
    clean: bool = True,
    bins: int = 0,
    prok: bool = False,
    filterseq: int = 1,
    strand: str = "both",
    order: int = 4,
    order_non: int = 2,
    gcode: str = "1",
    motifopt: str = "1",
    width: int = 12,
    prestart: int = 6,
    fixmotif: int = 1,
    offover: bool = True,
    par: Optional[str] = None,
    gibbs: int = 3,
    test: bool = False,
    identity: float = 0.99,
    maxitr: int = 10,
    verbose: int = 0,
    version: bool = False,
) -> None:
    if verbose == 1:
        logging.basicConfig(filename="pygmst.log", filemode="w", level=logging.WARN)
    elif verbose == 2:
        logging.basicConfig(filename="pygmst.log", filemode="w", level=logging.INFO)
    elif verbose == 3:
        logging.basicConfig(filename="pygmst.log", filemode="w", level=logging.DEBUG)

    motif = True if motifopt == "1" else False  # for compatibility
    fixmotif = True if fixmotif == 1 else False  # for compatibility

    # Not entirely sure this is the best way to include these in the module
    # but I am ignorant of how else to do so
    current_module = __import__(__name__)
    if current_module.__name__ != "__main__":
        working_path = os.path.dirname(sys.modules["pygmst"].__file__)
    else:
        working_path = f"{os.getcwd()}/pygmst"

    if output is None:
        base = os.path.basename(input)
        output = os.path.splitext(base)[0]
        if format == "LST":
            output = f"{output}.lst"
        elif format == "GFF":
            output = f"{output}.gff"
        if fnn:
            fnn_out = f"{output}.fnn"
        if faa:
            faa_out = f"{output}.faa"
    else:
        if fnn:
            fnn_out = f"{output}.fnn"
        if faa:
            faa_out = f"{output}.faa"
    if prok:
        bins = 1
        filterseq = 0
        order = 2
        order_non = 2
        offover = False
        gcode = "11"
        fixmotif = False
        prestart = 40
        width = 6
    if par is None:
        par = f"{working_path}/genemark/par_{gcode}.default"
    if version:
        print(f"{__version__}")
    if verbose:
        # TODO: add actual logging
        # log = logging.basicConfig(filename='gmst.log', level=logging.INFO)
        pass

    with tempfile.TemporaryDirectory() as tmpdir:

        # GeneMark.hmm gene finding program <gmhmmp>; version 2.14
        hmm = f"{working_path}/genemark/gmhmmp"

        # sequence parsing tool <probuild>; version 2.11c
        build = f"{working_path}/genemark/probuild"

        # gibbs sampler - from http://bayesweb.wadsworth.org/gibbs/gibbs.html ; version 3.10.001  Aug 12 2009
        # this doesn't appear to be available from that source?
        # possible replacement at https://github.com/Etschbeijer/GibbsSampler ?
        gibbs3 = f"{working_path}/genemark/Gibbs3"

        # use <probuild> with GeneMarkS parameter file <$par>
        build = f"{build} --par"

        # set options for <gmhmmp>

        # switch gene overlap off in GeneMark.hmm; for eukaryotic intron-less genomes
        if offover:
            hmm = f"{hmm} -p 0"

        # set strand to predict
        if strand == "direct":
            hmm = f"{hmm} -s d "
        elif strand == "reverse":
            hmm = f"{hmm} -s r "

        list_of_temp: List[str] = list()
        GC = 0

        # ------------------------------------------------
        ## tmp solution: get sequence size, get minimum sequence size from --par <file>
        ## compare, skip iterations if short

        # run(f"{build} --clean_join {seq} --seq {seqfile} --log {logfile}")

        # this should end up as something like
        # 'probuild --par par_1.default --clean_join sequence --seq test.fa
        run(f"{build} {par} --clean_join {seq} --seq {seqfile}".split())
        # list_of_temp.extend([seqfile])

        with open(seqfile, "r") as _:
            sequence_size = len(_.read())

        command = run(args=["grep", "MIN_SEQ_SIZE", par], capture_output=True)
        minimum_sequence_size = int(
            re.findall(
                pattern=r"\s*--MIN_SEQ_SIZE\s+(\d+)",
                string=str(command.stdout, "utf=8"),
            )[0]
        )

        do_iterations = 1

        if sequence_size < minimum_sequence_size:
            do_iterations = 0

        if do_iterations > 0:
            # ------------------------------------------------
            # clustering

            # initial prediction using MetaGeneMark model
            list_of_temp.extend([meta_out, f"{meta_out}.fna"])

            # get GC of whole sequence for each sequence"
            gc_out = f"{tmpdir}/{meta_out}.feature"

            # form of! 'probuild --par par_1.default --stat_fasta --seq test.fa > initial.meta.list.feature'
            gc_cmd = f"{build} {par} --stat_fasta --seq {seqfile}"
            with open(gc_out, "w") as gc_capture:
                logging.debug(gc_cmd)
                gc_capture.write(
                    str(run(gc_cmd.split(), capture_output=True).stdout, "utf-8")
                )

            list_of_temp.extend([gc_out])

            # determine bin number and range
            # TODO: I don't believe the cutoff values returned
            bin_num, cutoffs, seq_GC = cluster(
                feature_f=gc_out, clusters=bins, min_length=10000
            )

            # bin_num = 3 #debugging
            # Log("bin number = $bin_num\n");
            # Log( "GC range = ".join(",",@$cutoffs)."\n" );

            # ------------------------------------------------
            # training

            models: List[str] = list()
            seqs: List[str] = list()
            # my %handles; # file handles for n bins.
            if bin_num == 1:
                logging.info(
                    f"train(input_seq=seqfile, seq={seq}, motif={motif},fixmotif={fixmotif},order={order},order_non={order_non},start_prefix={start_prefix},gibbs_prefix={gibbs_prefix},prestart={prestart},width={width},build_cmd={build},hmm_cmd={hmm},par={par},maxitr={maxitr},identity={identity},gibbs3={gibbs3},workpath={working_path}"
                )
                final_model, new_tmp_files = train(
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
                    workpath=working_path,
                    tmpdir=tmpdir,
                )
                logging.info(
                    f"final_model: {final_model}, new_tmp_files: {new_tmp_files}"
                )
                list_of_temp.extend(new_tmp_files)
                # -----------------------------
                # make a GC file for the input file
                # open NEWINPUT, ">", $newseq;
                # -----------------------------
                # read input sequences
                try:
                    FA = pyfaidx.Fasta(seqfile)
                    seq_GC_entries = [_.lstrip(">") for _ in seq_GC.keys()]
                    for read in FA:
                        if read.long_name not in seq_GC_entries:
                            seq_GC[f">{read.long_name}"] = int(getGC(read[:].seq))
                except:
                    print(f"Cannot open {seqfile}")

            else:
                # create sequence file for each bin
                logging.info("Binning the input sequence")
                handles = SortedDict()
                for i in range(1, bin_num + 1):
                    with open(file=f"{tmpdir}/seq_bin_{i}", mode="w") as fh:
                        seqs.extend([f"{tmpdir}/seq_bin_{i}"])
                        handles[i] = f"{tmpdir}/seq_bin_{i}"
                        logging.info(f"Created {tmpdir}/seq_bin_{i}")

                # read input sequences
                try:
                    FA = pyfaidx.Fasta(seqfile)
                    seq_GC_entries = [_.lstrip(">") for _ in seq_GC.keys()]
                    for read in FA:
                        read_header = read.long_name
                        if read.long_name not in seq_GC_entries:
                            seq_GC[read_header] = int(getGC(read[:].seq))

                        # decide which bin the sequence belongs to
                        if bin_num == 2:
                            if seq_GC[read_header] <= cutoffs[0]:
                                bin = 1
                            else:
                                bin = 2
                        else:
                            if seq_GC[read_header] <= cutoffs[0]:
                                bin = 1
                            elif seq_GC[read_header] <= cutoffs[1]:
                                bin = 2
                            else:
                                bin = 3
                        # logging.info(f"Placed {read_header} into bin {bin}")

                        # output to corresponding output bin file
                        with open(file=handles[bin], mode="w") as fh:
                            fh.writelines(
                                f"{read_header}\t[gc={seq_GC[read_header]}\n{read[:].seq}]\n"
                            )
                except:
                    logging.info(f"Cannot open {seqfile}")
                    print(f"Cannot open {seqfile}")
                # train
                for i in range(bin_num):
                    current_model, new_tmp_files = train(
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
                        workpath=working_path,
                        bin_num=i,
                        tmpdir=tmpdir,
                    )
                    logging.info(f"train() returned: {current_model}")
                    if current_model:
                        models.extend([current_model])
                    list_of_temp.extend(new_tmp_files)
                # combine individual models to make the final model file
                final_model = combineModel(models, cutoffs)
                logging.info(f"combineModel() returned: {final_model}")

            run(f"cp {final_model} {out_name}".split())
            logging.info(f"Ran 'cp {final_model} {out_name}' but not sure why?")
            if out_name != final_model:
                list_of_temp.extend([final_model])

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

        list_of_temp.extend(["with_seq.out"])

        if do_iterations:
            if motif:
                command = (
                    f"{hmm} -r -m {out_name} -o {output} {format_gmhmmp} {seqfile}"
                )
                result = run(command.split(), stdout=PIPE, stderr=STDOUT)
                logging.debug(command)
                logging.debug(f"line: 533 {str(result.stdout, 'utf-8')}")
            else:
                # no moitf option specified
                command = f"{hmm} -m {out_name} -o {output} {format_gmhmmp} {seqfile}"
                result = run(command.split(), stdout=PIPE, stderr=STDOUT)
                logging.debug(command)
                logging.debug(f"line: 539 {str(result.stdout, 'utf-8')}")
        else:
            meta_model = resource_filename("pygmst", "genemark/MetaGeneMark_v1.mod")
            # no iterations - use heuristic only
            command = f"{hmm} -m {meta_model} -o {output} {format_gmhmmp} {seqfile}"
            result = run(command.split(), stdout=PIPE, stderr=STDOUT)
            logging.debug(command)
            logging.debug(f"line: 546 {str(result.stdout, 'utf-8')}")

        # this seems stupid dangerous.
        # TODO: replace with using tempfile.TemporaryDirectory
        # if clean:
        #     for _ in list_of_temp:
        #         run(["rm", "-f", _])
        # Final list was:
        # ['initial.meta.lst', 'initial.meta.lst.fna',
        #  '/tmp/tmpv1seq_gr/initial.meta.list.feature',
        #  'itr_0.lst', 'itr_1.mod', 'itr_1.lst', 'itr_2.mod',
        #  'itr_2.lst', 'itr_3.mod', 'itr_3.lst', 'itr_3.mod',
        #  'with_seq.out']
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
    workpath: str,
    bin_num: int = 0,
    tmpdir: Optional[str] = None,
) -> Tuple[Optional[str], List[str]]:
    tmp_files: List[str] = list()
    logging.info("Beginning training")
    # ------------------------------------------------
    # prepare sequence
    # build_cmd should have the form of! `probuild --par par_1.default`
    # so this makes the below
    # `probuild --par par_1.default --clean_join sequence --seq test.fa`
    # which seems kind of redundant, but what do I know

    # TODO: we need to check that the input_sequence isn't empty here
    command = f"{build_cmd} {par} --clean_join {seq} --seq {input_seq}"
    result = run(command.split(), stdout=PIPE, stderr=STDOUT)
    logging.debug(command)
    logging.debug(f"line: 596 {str(result.stdout, 'utf-8')}")
    # TODO: another thing to restore when we figure out logging.
    # run([build_cmd, "--clean_join", seq, "--seq", input_seq, "--log", logfile])

    # ------------------------------------------------
    # tmp solution: get sequence size, get minimum sequence size from --par <file>
    # compare, skip iterations if short

    with open(seq, "r") as _:
        sequence_size = len(_.read())

    min_size_find_cmd = run(
        args=f"grep MIN_SEQ_SIZE {par}".split(), capture_output=True
    )
    minimum_sequence_size = int(
        re.findall(
            pattern=r"\s*--MIN_SEQ_SIZE\s+(\d+)",
            string=str(min_size_find_cmd.stdout, "utf=8"),
        )[0]
    )

    do_iterations = True
    iterated = bool()
    itr = 0

    logging.info(f"minimum_sequence_size: {minimum_sequence_size}")
    logging.info(f"sequence_size: {sequence_size}")
    if sequence_size < minimum_sequence_size:
        # form of! `probuild --par par_1.default --clean_join sequence --seq test.fa --MIN_CONTIG_SIZE 0 --GAP_FILLER
        command = f"{build_cmd} {par} --clean_join {tmpdir}/{seq} --seq {input_seq} --MIN_CONTIG_SIZE 0 --GAP_FILLER"
        result = run(command.split(), stdout=PIPE, stderr=STDOUT)
        logging.debug(command)
        logging.debug(f"line: 628 {str(result.stdout, 'utf-8')}")

        init_mod = resource_filename("pygmst", "genemark/MetaGeneMark_v1.mod")
        next_item = f"{tmpdir}/{hmmout_prefix}{itr}_bin_{bin_num}{hmmout_suffix}"
        command = f"{hmm_cmd} -m {init_mod} -o {next_item} {tmpdir}/{seq}"
        result = run(command.split(), stdout=PIPE, stderr=STDOUT)
        logging.debug(command)
        logging.debug(f"line: 635 {str(result.stdout, 'utf-8')}")

        mod = f"{tmpdir}/bin_{bin_num}{mod_suffix}"
        command = f"{build_cmd} {par} --mkmod {mod} --seq {tmpdir}/{seq} --geneset {next_item} --ORDM {order} --order_non {order_non} --revcomp_non 1"

        if motif and not fixmotif:
            command = f"{command} --pre_start {start_seq} --PRE_START_WIDTH {prestart}"
        elif motif and fixmotif:
            command = (
                f"{command} --fixmotif --PRE_START_WIDTH {prestart} --width {width}"
            )

        results = run(command.split(), stdout=PIPE, stderr=STDOUT)
        logging.debug(command)
        logging.debug(f"line: 649 {str(result.stdout, 'utf-8')}")
        return mod, tmp_files

    # Log("do_iterations = $do_iterations\n")

    # ------------------------------------------------
    # run initial prediction
    next_item = f"{hmmout_prefix}{itr}_bin_{bin_num}{hmmout_suffix}"
    print("run initial prediction")

    # form of! `gmhmmp -s d sequence -m MetaGeneMark_v1.mod -o itr_{itr}.lst`
    mod = resource_filename("pygmst", "genemark/MetaGeneMark_v1.mod")
    command = f"{hmm_cmd} -m {mod} -o {next_item} {seq}"
    result = run(command.split(), stdout=PIPE, stderr=STDOUT)
    logging.debug(command)
    logging.debug(f"line: 664 {str(result.stdout, 'utf-8')}")

    # TODO: tempfile
    tmp_files.extend([next_item])

    # enter iterations loop
    if do_iterations:
        logging.info("entering training iteration loop")
    while do_iterations:
        iterated = True
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
        results = run(command.split(), stdout=PIPE, stderr=STDOUT)
        logging.debug(command)
        logging.debug(f"line: 695 {str(result.stdout, 'utf-8')}")

        tmp_files.extend([mod])

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
            result = run(command.split(), stdout=PIPE, stderr=STDOUT)
            logging.debug(command)
            logging.debug(f"line: 711 {str(result.stdout, 'utf-8')}")
            tmp_files.extend([start_seq])

            logging.info("make prestart model")
            # TODO: logfile
            # run([build_cmd, "--gibbs", gibbs_out, "--mod", mod, "--seq", start_seq, "--log", logfile])
            # form of! `probuild --par par_1.default --gibs gibbs_out.{itr} --mod itr_{itr}.mod --seq startseq.{itr}
            command = (
                f"{build_cmd} {par} --gibbs {gibbs_out} --mod {mod} --seq {start_seq}"
            )
            result = run(command.split(), stdout=PIPE, stderr=STDOUT)
            tmp_files.extend([gibbs_out])
            logging.debug(command)
            logging.debug(f"line: 724 {str(result.stdout, 'utf-8')}")

        prev = next_item
        # next_item = itr_{itr}.lst
        next_item = f"{hmmout_prefix}{itr}_bin_{bin_num}{hmmout_suffix}"

        # form of! `gmhmmp -s d sequence -m itr_{itr}.mod -o itr_{itr}.lst
        command = f"{hmm_cmd} {seq} -m {mod} -o {next_item}"

        # form of! `gmhmmp -s d sequence -m itr_{itr}.mod -o itr_{itr}.lst -r
        if motif:
            command += " -r"

        logging.info(f"prediction, iteration: {itr}")
        result = run(command.split(), stdout=PIPE, stderr=STDOUT)
        logging.debug(command)
        logging.debug(f"line: 740 {str(result.stdout, 'utf-8')}")
        tmp_files.extend([next_item])

        # `probuild --par par_1.default --compare --source itr_{itr}.lst --target itr_{itr-1}.lst
        command = f"{build_cmd} {par} --compare --source {next_item} --target {prev}"
        # &Log( "compare:\n" . $command . "\n" );

        results = run(command.split(), stdout=PIPE, stderr=STDOUT)
        logging.debug(command)
        logging.debug(f"line: 749 {str(result.stdout, 'utf-8')}")

        diff = str(results.stdout, "utf-8").strip("\n")
        logging.info(f"iteration {itr}, bin {bin_num} difference: {diff}")
        # &Log( "compare $prev and $next_item: $diff\n" );

        if float(diff) >= identity:
            logging.info(f"Stopped iterations on identity: {diff}")
            do_iterations = False
        if itr == maxitr:
            logging.info(f"Stopped iterations on maximum number: {maxitr}")
            return mod, tmp_files
    logging.info(f"return value of train: {mod}")
    return mod, tmp_files


def cluster(
    feature_f: str, clusters: int, min_length: int = 10000
) -> Tuple[int, List[int], Dict[str, int]]:  # $gc_out, $bins
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
            logging.debug(f"line 781: text.group(0) = {text.group(0)}")
            if text:
                header = text.group(1)  # Reference name
                logging.debug(f"line 784: header = {header}")
                length = int(text.group(2))  # length of sequence?
                logging.debug(f"line 786: header = {length}")
                seqGC = int(text.group(3))  # must be GC percentage
                logging.debug(f"line 788: header = {seqGC}")

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
        if previous < total_length / 3 * 2 and gc_hash[key] >= total_length / 3 * 2:
            two_third = key
        if previous < total_length / 2 and gc_hash[key] >= total_length / 2:
            one_half = key
        previous = gc_hash[key]
        # TODO: uncomment when we have logging fixed
        # log.info(f"({one_third})->({gc_hash[one_third]})\n")
        # log.info(f"({one_half})->({gc_hash[one_half]})\n")
        # log.info(f"({two_third})->({gc_hash[two_third]})\n")

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


# -----------------------------------------------
# concatenate model files into one in MGM format:
# starts with "__GC" and ends with "end"
# -----------------------------------------------


def combineModel(mod: List[str], cut_offs: List[int], minGC: int = 30, maxGC: int = 70):
    # change the min and max GC value of cut_offs to the minGC and maxGC used by gmhmmp
    cut_offs[0] = minGC
    cut_offs[-1] = maxGC

    b = 1
    data = str()
    with open(file="final_model", mode="w") as model:
        for i in range(minGC, maxGC + 1):
            model.write(f"__GC{i}\t\n")
            if i == cut_offs[b] or i == maxGC:
                with open(file=mod[b - 1], mode="r") as fh:
                    for line in fh:
                        if "NAME" in line:
                            line = line.strip("\n")
                            line = f"{line}_GC<={i}\n"
                        data += line
                model.write(data)
                model.write("end \t\n\n")
                if b > len(mod):
                    break
                b += 1
    return "final_model"


if __name__ == "__main__":
    if "--verbose" in sys.argv:
        verbose = True
    sys.exit(main())  # pragma: no cover
