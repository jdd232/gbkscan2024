#!/usr/bin/perl
use 5.14.0;
use utf8;
use Carp;
use lib qw(/home/sco/perllib);
use File::Basename;
use Getopt::Long;
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
  tablistV tablistVE linelistV linelistVE tablistH linelistH
  tablistER tablistVER linelistER linelistVER tabhashER tabhashVER csvsplit);
use File::Spec;
use File::Temp qw(tempfile tempdir);
use DBI;
use Cwd;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
my $command    = join( " ", $0, @ARGV );
my $workingdir = getcwd();

# {{{ Getopt::Long
my $outdir;
my $indir;
my $fofn;
my $outex; # extension for the output filename when it is derived on infilename.
my $errfile;
my $allowed_dist = 200;    # max distance from gene.
my $gbkfn;
my $allowed = 2;             # mismatches allowed.
my $left    = q/taaa/;
my $right   = q/gccgataa/;
my $gaplen  = 16;
my $runfile;
my $beginM = 1;
my $protfna_wanted = 1;
my $outfile;
my $testCnt = 0;
our $verbose;
my $skip = 0;
my $help;
GetOptions(
  "outfile:s"            => \$outfile,
  "gbkfn|gbkout=s"       => \$gbkfn,
  "protfna!"             => \$protfna_wanted,
  "beginmet!"            => \$beginM,
  "allowed|mismatches:i" => \$allowed,
  "maxdistance:i"        => \$allowed_dist,
  "left:s"               => \$left,
  "right:s"              => \$right,
  "dirout:s"             => \$outdir,
  "gaplen:i"             => \$gaplen,
  "indir:s"              => \$indir,
  "fofn:s"               => \$fofn,
  "extension:s"          => \$outex,
  "errfile:s"            => \$errfile,
  "runfile:s"            => \$runfile,
  "testcnt:i"            => \$testCnt,
  "skip:i"               => \$skip,
  "verbose"              => \$verbose,
  "help"                 => \$help
);

# }}}

# {{{ POD

=head1 Name

gbk_scan.pl

=head2 Examples

 gbk_scan.pl -gbkfn out/RradSPS_search.gbk \
 -left taaa -right gccgataa -gaplen 16 -maxdistance 200 \
 -outfile out/RradSPS_search.txt \
 -- gbk/GCA_000661895.1_ASM66189v1_genomic.gbff

 print_header.pl -outfile $outdir/RradSPS_search.txt -- mark.header


=head2 Description

Search -10 and -35 with mismatches allowed on both side.

=head2 Notes

=head3 Matches

If a match falls within a gene then it is not reported (sub lir)

If a match is more than I<$allowed_dist> from the CDS start in the right
direction then it is not reported (sub lir).

=head3 Output files

Tabular output is written to the file specified by the C<-outfile>
option.

If a genbank output file is specified by the C<-gbkfn> or C<-gbkout>
option then as many genbank files are written as there are entries in
the input genbank file. In these genbank files the locations matches are
tagged a I<protein_binding> with an attribute added so that they will be
rendered in red by Artemis.

For example, if you specified S<C<-gbkfn RradSPS_search.gbk>> and the
input genbank file has 4 entries in it, then you will get the
following four genbank files

* RradSPS_search.gbk
* RradSPS_search_p1.gbk
* RradSPS_search_p2.gbk
* RradSPS_search_p3.gbk

This is to allow you to view the different entries separately in Artemis.

=cut

# }}}

if ($help) {
  exec("perldoc $0");
  exit;
}

# {{{ open the errfile
if ($errfile) {
  open( ERRH, ">", $errfile );

  # print(ERRH "$0", "\n");
  close(STDERR);
  open( STDERR, ">&ERRH" );
}

# }}}

# {{{ Temporary directory and template.
# my($tmpfh, $tmpfn)=tempfile($template, DIR => $tempdir, SUFFIX => '.tmp');
# somewhere later you need to do this
# unlink($tmpfn);
# unlink(glob("$tmpfn*"));
# }}}

# {{{ populate @infiles
my @infiles;
if ( -e $fofn and -s $fofn ) {
  open( FH, "<", $fofn );
  while ( my $line = readline(FH) ) {
    chomp($line);
    if ( $line =~ m/^\s*\#/ or $line =~ m/^\s*$/ ) { next; }
    my $fn;
    if ($indir) {
      $fn = File::Spec->catfile( $indir, $line );
    }
    else {
      $fn = $line;
    }

    push( @infiles, $fn );
  }
  close(FH);
}
else {
  @infiles = @ARGV;
}

# }}}

# {{{ outdir and outfile business.
my $ofh;
my $idofn = 0;    # Flag for input filename derived output filenames.
if ($outfile) {
  my $ofn;
  if ($outdir) {
    unless ( -d $outdir ) {
      unless ( mkdir($outdir) ) {
        croak("Failed to make $outdir. Exiting.");
      }
    }
    $ofn = File::Spec->catfile( $outdir, $outfile );
  }
  else {
    $ofn = $outfile;
  }
  open( $ofh, ">", $ofn );
}
else {
  open( $ofh, ">&STDOUT" );
}
select($ofh);
# }}}

my @header = qw(serial  accession  location  sequence
  mismatches  strand  dist2gene locus_tag  protein_id  product);
$header[0] = "#" . $header[0];
tablist(@header);


my $tempdir  = qw(/tmp);
my $template = "featsq3XXXXX";
my ( $sq3fh, $sq3fn ) =
  tempfile( $template, DIR => $tempdir, SUFFIX => '.sqlite3' );
close($sq3fh);
my $handle = DBI->connect( "DBI:SQLite:dbname=$sq3fn", '', '' );

my $infile = shift(@infiles);

$left  = lc($left);
$right = lc($right);

my @left  = split( //, $left );
my @right = split( //, $right );

my $ifh;
if ( $infile =~ m/\.gz$/ ) {
  $ifh = gzfh($infile);
}
else {
  open( $ifh, "<$infile" ) or croak("Could not open $infile");
}
my $seqio = Bio::SeqIO->new( -fh => $ifh );

my $protfnafh; my $protfnaout;
my $protfaafh; my $protfaaout;
if ($gbkfn) {
  my ( $noex, $dir, $ext ) = fileparse( $gbkfn, qr/\.[^.]*/ );
  my $protfnafn = File::Spec->catfile( $dir, $noex . "_prot.fna" );
  open($protfnafh, ">", $protfnafn );
  $protfnaout = Bio::SeqIO->new( -fh => $protfnafh, -format => 'fasta' );
  my $protfaafn = File::Spec->catfile( $dir, $noex . "_prot.faa" );
  open($protfaafh, ">", $protfaafn );
  $protfaaout = Bio::SeqIO->new( -fh => $protfaafh, -format => 'fasta' );
}

my $entry = 0;
while ( my $seqobj = $seqio->next_seq() ) {
  my $seqlen = $seqobj->length();
  my $lastdx = $seqlen - 1;
  my $revobj = $seqobj->revcom();

  $handle->begin_work();

  # Make a fresh sql table for each entry in the gbk file.
  my $dropstr = qq/drop table if exists featemp/;
  unless ( $handle->do($dropstr) ) {
    croak("Failed: $dropstr");
  }

  # {{{ Create temporary table featemp.
  my $creatstr = qq/create temporary table featemp (/;
  $creatstr .= qq/ lt text not null,/;
  $creatstr .= qq/ protid text not null,/; # If there is no protid we will put 'no_protein_id'.
  $creatstr .= qq/ fstart integer not null,/;
  $creatstr .= qq/ fend integer not null,/;
  $creatstr .= qq/ fstrand integer not null,/;
  $creatstr .= qq/ product text/;
  $creatstr .= qq/)/;

  unless ( $handle->do($creatstr) ) {
    linelistE($creatstr);
    croak("Failed to create table feat");
  }

  # }}}

  # {{{ Temporary table is populated here.
  for my $feat ( $seqobj->all_SeqFeatures() ) {
    if ( $feat->primary_tag() eq 'CDS' ) {
      my ( $lt, $protid, $start, $end, $strand, $product );
      if ( $feat->has_tag("locus_tag") ) {
        my @temp = $feat->get_tag_values("locus_tag");
        $lt     = $temp[0];
        $start  = $feat->start();
        $end    = $feat->end();
        $strand = $feat->strand();
        if ( $feat->has_tag("product") ) {
          my @temp = $feat->get_tag_values("product");
          $product = join( " ", @temp );
          $product =~ s/'//g;
        }
        if ( $feat->has_tag("protein_id") ) {
          my @temp = $feat->get_tag_values("protein_id");
          $protid = $temp[0];
        } else { $protid = "no_protein_id"; }
        my @invals = (
          $handle->quote($lt), $handle->quote($protid),
          $start, $end, $strand, $handle->quote($product)
        );
        my $instr = qq/insert into featemp values(/;
        $instr .= join( ",", @invals );
        $instr .= qq/)/;
        if ( $handle->do($instr) ) {
          my $noop = 1;
        }
        else {
          linelistE($instr);
        }
      }
    }
  }

  # }}}

  my ( $gbfh, $seqout );
  if ($gbkfn) {
    my ( $noex, $dir, $ext ) = fileparse( $gbkfn, qr/\.[^.]*/ );
    if ($entry) { $noex .= "_p$entry"; }
    my $gbfn = File::Spec->catfile( $dir, $noex . $ext );
    open( $gbfh, ">", $gbfn );
    $seqout = Bio::SeqIO->new( -fh => $gbfh, -format => 'genbank' );
  }
  my $id     = $seqobj->display_id();
  my $seq    = $seqobj->seq();
  my $revseq = $revobj->seq();
  my $ssdx   = 0;
  my @strand = ( 1, -1 );
  my @seqobj = ( $seqobj, $revobj );
  my @seq    = ( $seq, $revseq );
  my @finds;
  linelistE("Scanning $id");

  for my $strandseq (@seq) {
    my $strand    = $strand[$ssdx];
    my $strandobj = $seqobj[$ssdx];
    my $gps       = 0;
    my $subseq;
    while ( $subseq =
      substr( $strandseq, $gps, length($left) + $gaplen + length($right) ) )
    {
      if ( length($subseq) < length($left) + $gaplen + length($right) ) {
        last;
      }
      my $lside  = substr( $subseq, 0, length($left) );
      my $spacer = substr( $subseq, length($left), $gaplen );
      my $rside  = substr( $subseq, -( length($right) ) );
      my ( $lmis, $rmis, $miscnt ) = mismatches( $lside, $rside );
      if ( $miscnt <= $allowed ) {
        push( @finds, [ $gps, $strand, $subseq, "$lmis:$rmis" ] );

        # tablistE($id, $gps, $strand, "$lside:$lmis",
        # $spacer, "$rside:$rmis", $miscnt);
      }
      $gps += 1;

      # if(scalar(@finds) >= 3) { last; }
      unless ( $gps % 300_000 ) {
        linelistER("$gps           of $seqlen");
      }
    }
    $ssdx += 1;
  }

  my @ffinds;
  for my $lr (@finds) {
    my $gps    = $lr->[0];
    my $strand = $lr->[1];
    my $subseq = $lr->[2];
    my $miss   = $lr->[3];

    # unless ($miss =~ m/^0:/) { next; } ### testing only.
    if ( $strand == -1 ) {
      $gps += ( length($subseq) - 1 );
      $gps = $lastdx - $gps;
      push( @ffinds, [ $gps, $gps + ( length($subseq) - 1 ), $strand, $miss ] );
    }
    else {
      push( @ffinds, [ $gps, $gps + ( length($subseq) - 1 ), $strand, $miss ] );
    }
  }
  tablistE( "\nCount ignoring context in $id:", scalar(@ffinds), "\n" );
  my $outserial = 1;
  my %donelt;
  for my $lr ( sort sorter @ffinds ) {
    my $start  = $lr->[0] + 1;
    my $end    = $lr->[1] + 1;
    my $strand = $lr->[2];
    my $miss   = $lr->[3];
    my ( $lt, $protid, $product, $dist ) = lir(
      start  => $start,
      end    => $end,
      strand => $strand
    );
    if ( $lt and $dist <= $allowed_dist ) {
      my $withFlank = with_flank( $seqobj, $start, $end, $strand );

  # tablistE($outserial, $id, $fpos + 1, $found, $strand, $dist, $lt, $product);
      tablist(
        $outserial, $id,   $start, $withFlank, $miss,
        $strand,    $dist, $lt, $protid, $product
      );
      $outserial += 1;
      if ($gbkfn) {
        my $tagnote = "$left N$gaplen $right $dist nt from $lt.";
        my $whifeat = Bio::SeqFeature::Generic->new(
          -primary => 'protein_bind',
          -start   => $start,
          -end     => $end,
          -strand  => $strand,
          -tag     => {
            'note'   => $tagnote,
            'colour' => "255 0 0"
          }
        );
        $seqobj->add_SeqFeature($whifeat);
      }
      if($protfna_wanted) {
        unless(exists($donelt{$lt})) {
          my %cdo = cdsobjects(seqobj => $seqobj, lt => $lt);
          my $cdsfnaobj = $cdo{nt};
          my $cdsfaaobj = $cdo{prot};
          $protfnaout->write_seq($cdsfnaobj);
          $protfaaout->write_seq($cdsfaaobj);
          $donelt{$lt} = 1;
        }
      }
    }
  }
  if ($gbkfn) {
    $seqout->write_seq($seqobj);
    close($gbfh);
  }
  $entry += 1;
  $handle->commit();
}
close($ifh);
close($protfnafh);
close($protfaafh);
exit;

sub sorter {
  return ( $a->[0] <=> $b->[0] );
}

# {{{ sub cdsobjects. Returns a hash with keys prot and nt.
sub cdsobjects {
  my %args = @_;
  my $lt = $args{lt};
  my $seqobj = $args{seqobj};
  my $qstr = qq/select fstart, fend, fstrand, product from featemp where lt = '$lt'/;
  my ($start, $end, $strand, $product) = $handle->selectrow_array($qstr);
  my $protobj;
  my $ntobj;
  if($strand == -1) {
    $ntobj = $seqobj->trunc($start, $end)->revcom();
  }
  else {
    $ntobj = $seqobj->trunc($start, $end);
  }
  $ntobj->display_id($lt);
  $ntobj->description($product);
  $protobj = $ntobj->translate();
  if($beginM) {
    my $aaseq = $protobj->seq();
    $aaseq =~ s/^./M/;
    $protobj->seq($aaseq);
  }
  return(prot => $protobj, nt => $ntobj);      
}
# }}}


# {{{ sub revcom
sub revcom {
  my $seq    = shift(@_);
  my $obj    = Bio::Seq->new( -seq => $seq );
  my $revobj = $obj->revcom();
  return ( $revobj->seq() );
}

# }}}

# {{{ sub mismatches
sub mismatches {
  my $lside = lc( shift(@_) );
  my $rside = lc( shift(@_) );
  my $lmis  = 0;
  my $rmis  = 0;
  my @lside = split( //, $lside );
  my @rside = split( //, $rside );
  while ( my ( $dx, $nt ) = each(@lside) ) {
    if ( $left[$dx] ne $nt ) { $lmis += 1; }
  }
  while ( my ( $dx, $nt ) = each(@rside) ) {
    if ( $right[$dx] ne $nt ) { $rmis += 1; }
  }

  my $miscnt = $lmis + $rmis;
  return ( $lmis, $rmis, $miscnt );
}

# }}}

# {{{ sub with_flank
# There is an unacceptable level of hard-coding in this
# subroutine. Some improvement is possible but it will
# remain messy. So I am leaving it as it is.
sub with_flank {
  my $seqobj = shift(@_);
  my $start  = shift(@_);
  my $end    = shift(@_);
  my $strand = shift(@_);
  my $retstr;
  if ( $strand > 0 ) {
    $start -= 10;
    $end = $start + 56;
    my $section = $seqobj->subseq( $start, $end );
    $retstr = $section;
  }
  else {
    $end   = $end + 10;
    $start = $end - 56;
    my $section = $seqobj->subseq( $start, $end );
    my $tempobj = Bio::Seq->new( -seq => $section );
    my $revobj  = $tempobj->revcom();
    my $revseq  = $revobj->seq();
    $retstr = $revseq;
  }
  substr( $retstr, 10 + length($left),          0 ) = ":";
  substr( $retstr, 10 + length($left) + 17,     0 ) = ":";
  substr( $retstr, 10 + length($left) + 17 + 9, 0 ) = ":";
  return ($retstr);
}

# }}}

# {{{ sub lir
sub lir {
  my %argv   = @_;
  my $start  = $argv{start};
  my $end    = $argv{end};
  my $strand = $argv{strand};
  my $qstr1  = qq/select lt from featemp where fstart <= $start/;
  $qstr1 .= qq/ and fend >= $end/;
  my ($within) = $handle->selectrow_array($qstr1);
  if ($within) {
    return ();
  }

  if ( $strand == -1 ) {
    my $qleft =
    qq/select lt, protid, fstart, fend, fstrand, product from featemp where fend = /;
    $qleft .= qq/(select max(fend) from featemp where fend < $start/;
    $qleft .= qq/ and fstrand = $strand)/;
    my ( $lt, $protid, $fstart, $fend, $fstrand, $product ) =
      $handle->selectrow_array($qleft);
    my $dist = ( $start - $fend ) - 1;
    return ( $lt, $protid, $product, $dist );
  }

  elsif ( $strand == 1 ) {
    my $qright =
    qq/select lt, protid, fstart, fend, fstrand, product from featemp where fstart = /;
    $qright .= qq/(select min(fstart) from featemp where fstart > $end/;
    $qright .= qq/ and fstrand = $strand)/;
    my ( $lt, $protid, $fstart, $fend, $fstrand, $product ) =
      $handle->selectrow_array($qright);
    my $dist = ( $fstart - $end ) - 1;
    return ( $lt, $protid, $product, $dist );
  }

}

# }}}

# {{{ sub gzfh.
sub gzfh {
  my $ifn = shift(@_);
  my $ifh;
  $ifh = tempfile();
  if ( gunzip $infile => $ifh ) {
    seek( $ifh, 0, 0 );
  }
  else {
    close($ifh);
    linelistE("gunzip failed: $GunzipError");
  }
  return ($ifh);
}

# }}}

# Multiple END blocks run in reverse order of definition.
END {
  close(STDERR);
  close(ERRH);
  close($ofh);
  $handle->commit();
  $handle->disconnect();
  unlink($sq3fn);
}

__END__


