package Sco::Common;
use 5.14.0;
use Exporter qw(import);
use File::Temp qw(tempfile tempdir);
our @EXPORT_OK = qw(
tablist linelist tablistE linelistE tabhash tabhashE tabvals
tablistV tablistVE linelistV linelistVE tablistH linelistH
tablistER tablistVER linelistER linelistVER tabhashER tabhashVER
ymd ddmyhms csvsplit file2hash linelistserial csvlist inlist
gzfh
);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
# use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
#                    tablistV tablistVE linelistV linelistVE);

my @days = qw(Sun Mon Tue Wed Thu Fri Sat);

# {{{ subroutines tablist, linelist, tabhash and their *E versions.
# The E versions are for printing to STDERR.
# tabvals is also here.

sub linelistserial {
  my @in = @_;
  my $serial = 0;
  for my $el (@in) {
    $serial += 1;
    print("$serial\t$el\n");
  }
}

sub tablistVER {
  if($main::verbose) {
  my @in = @_;
  print(STDERR join("\t", @in), "\r");
  }
}


sub tablistVE {
  if($main::verbose) {
  my @in = @_;
  print(STDERR join("\t", @in), "\n");
  }
}

sub tablistV {
  if($main::verbose) {
  my @in = @_;
  print(join("\t", @in), "\n");
  }
}

sub tablist {
  my @in = @_;
  print(join("\t", @in), "\n");
}

sub tablistER {
  my @in = @_;
  print(STDERR join("\t", @in), "\r");
}


sub tablistE {
  my @in = @_;
  print(STDERR join("\t", @in), "\n");
}


sub linelistV {
  if($main::verbose) {
  my @in = @_;
  print(join("\n", @in), "\n");
  }
}

sub linelistVE {
  if($main::verbose) {
  my @in = @_;
  print(STDERR join("\n", @in), "\n");
  }
}

sub linelistVER {
  if($main::verbose) {
  my @in = @_;
  print(STDERR join("\n", @in), "\r");
  }
}


sub linelist {
  my @in = @_;
  print(join("\n", @in), "\n");
}

sub linelistE {
  my @in = @_;
  print(STDERR join("\n", @in), "\n");
}

sub linelistER {
  my @in = @_;
  print(STDERR join("\n", @in), "\r");
}


sub tabhash {
  my %in = @_;
  for my $key (sort keys(%in)) {
    print(join("\t", $key, $in{$key}), "\n");
  }
}

sub tabhashE {
  my %in = @_;
  for my $key (sort keys(%in)) {
    print(STDERR join("\t", $key, $in{$key}), "\n");
  }
}

sub tabhashER {
  my %in = @_;
  for my $key (sort keys(%in)) {
    print(STDERR join("\t", $key, $in{$key}), "\r");
  }
}

sub tabvals {  # takes a hashref and a listref.
  my $hr = shift(@_);
  my $lr = shift(@_);
  my @temp;
  if(ref($lr)) {
    for my $key (@{$lr}) {
      push(@temp, $hr->{$key});
    }
  }
  else {
    for my $key (@{$hr}) {
      push(@temp, $hr->{$key});
    }
  }
  print(join("\t", @temp), "\n");
}

# }}}

# {{{ tablistH and linelistH
sub tablistH {
my $fh = shift(@_);
my @in = @_;
print($fh join("\t", @in), "\n");
}

sub linelistH {
  my $fh = shift(@_);
  my @in = @_;
  print($fh join("\n", @in), "\n");
}
# }}}

# {{{ sub ymd (year month date)
sub ymd {
my $time = shift(@_);
unless($time) { $time = time(); }
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);
$year += 1900;
my $month = $mon + 1;
return($year, $month, $mday);
}
# }}}

# {{{ sub ddmyhms (weekday date month year hour min sec)
sub ddmyhms {
my $time = shift(@_);
unless($time) { $time = time(); }
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);
$year += 1900;
my $xday = $days[$wday];
my $month = sprintf("%2d", $mon + 1);

my $rsec = sprintf("%2d", $sec);
my $rmin = sprintf("%2d", $min);
my $rhour = sprintf("%2d", $hour);

return($xday, $mday, $month, $year, $rhour, $rmin, $rsec);
}

# }}}

# {{{ sub csvsplit

=head2 Sub csvsplit

The string argument, typically a line read from a file,
is split by commas. Any commas inside double quotes are
ignored for the purpose of splitting.

This implementation is probably not very efficient. It is
certainly not flexible.


=cut

sub csvsplit {
  my $in = shift(@_);
  chomp($in);
  my $tok;
  my $qfl = 0;
  my @tok;
  for my $ch (split('', $in)) {
    if($ch eq '"') { $qfl = $qfl == 0 ? 1 : 0; }
    elsif($ch eq ",") {
      if($qfl) {
        $tok .= $ch;
      }
      else {
        push(@tok, $tok);
        $tok = '';
      }
    }
    else { $tok .= $ch; }
  }
  if($tok) { push(@tok, $tok); }
return(@tok);
}
# }}}

# {{{ sub csvlist
sub csvlist {
  my @in = @_;
  my @out;
  for my $in (@in) {
    if($in =~ m/[^+\-.0-9e]/i) {
      push(@out, qq("$in"));
    }
    else {
      push(@out, $in);
    }
  }
  print(join(",", @out), "\n");
}
# }}}

# {{{ sub file2hash

=head2 Sub file2hash

Takes a two column file (tab separated) and returns
a hash.

=cut

sub file2hash {
  my $ifn = shift(@_);
  open(my $ifh, "<", $ifn);
  my %rethash;
    while(my $line = readline($ifh)) {
      chomp($line);
      if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
      my @ll = split(/\t/, $line);
      $rethash{$ll[0]} = $ll[1];
    }
    close($ifh);
  return(%rethash);
}
# }}}

# {{{ sub inlist
sub inlist {
  my %in = @_;
  my $lr = $in{lr};
  my $el = $in{el};
  my $numeric = $in{numeric};
  if($numeric) {
    if(grep {$el == $_} @{$lr}) {
      return(1)
    }
  }
  else {
    if(grep {$el eq $_} @{$lr}) {
      return(1)
    }
  }
  return(0);
}
# }}}

# {{{ sub gzfh
sub gzfh {
  my $self = shift(@_);
  my $ifn = shift(@_);
  my $ifh = tempfile();
  if(gunzip $ifn => $ifh) {
    seek($ifh,0,0);
  }
  else {
    close($ifh);
    linelistE("gunzip failed: $GunzipError");
    return();
  }
  return($ifh);
}
# }}}


return(1);

