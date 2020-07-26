# In the event they are lost, this will generate the expected result for
# &combineModel()

use strict;
use Getopt::Long;
use FindBin qw($RealBin);
use File::Spec;

&combineModel();

sub combineModel{
  print \@_;
    my @mod = ("test_model_1.mod", "test_model_2.mod");
    my $mod = \@mod;
    my @cut_offs = {30, 70};
    my $minGC = 30;
    my $maxGC = 70;
#	my ($mod, $cut_offs) = @_;
    
    #change the min and max GC value of cut_offs to the minGC and maxGC used by gmhmmp
    
    my $b = 1;
    open MODEL, ">perl_final_model";
    for(my $i = $minGC; $i <=$maxGC; $i ++){
        print MODEL "__GC".$i."\t\n";
        if($i == $cut_offs[$b] || $i == $maxGC){
            open my $fh, '<', $mod->[$b-1] or die;
            #my $data = do { local $/; <$fh> };
            #close $fh;
            my $data;
            while(my $line = <$fh>){
                if($line =~ /NAME/){
                    chomp $line;
                    $line .= "_GC<=$i\n";
                }
                $data .= $line;
            }
            close $fh;
            print MODEL $data;
            print MODEL "end \t\n\n";
            last if( $b > scalar(@$mod));
            $b ++;
        }
    }
    close MODEL;
    return "final_model";
}