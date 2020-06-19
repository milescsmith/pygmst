my $handles;
my $bin_num = 3;

for(my $i = 1; $i <= $bin_num; ++$i){
    my $fh;
    open($fh, ">seq_bin_$i");
    push(@seqs, "seq_bin_$i");
    print($fh);
    print(@seqs);
    print("\n");
    $handles{$i} = $fh;
}