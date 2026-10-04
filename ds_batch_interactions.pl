#! /usr/bin/perl -w
#
# ds_batch_interactions.pl   (v2)
#
# Open this file from  File > Open  in Discovery Studio Visualizer.
# For EVERY *.pdb in  ~/Documents/ds_automation/pdb  it:
#   - sets a white background / publication quality,
#   - hides the protein ribbon (cartoon),
#   - finds the ligand automatically (any non-protein, non-water residue: LIG, UNL, *, 5T2 ...),
#   - creates the ligand non-bond interactions (H-bond, alkyl, pi-alkyl, pi-cation, pi-sulfur, metal ...),
#   - keeps ONLY the ligand + the residues that really interact with it (+ metals such as Zn),
#   - draws the ligand as ball & stick (element colours, custom carbon colour), residues as thin lines,
#   - labels each interacting residue (black text),
#   - saves   <name>_3D.png  and   <name>_2D.svg   in  ~/Documents/ds_automation/out
# A step-by-step report is written to  out/batch_log.txt
#
# Afterwards run  finish_images.py  (in a terminal) -> final/ folder with 600-dpi 2D + cropped 3D images.
#
use strict;
use MdmDiscoveryScript;
use DSCommands;

# =========================== SETTINGS (edit here) ===========================
my $HOME   = $ENV{HOME} || "/home/yaseen";
my $BASE   = "$HOME/Documents/ds_automation";
my $PDBDIR = "$BASE/pdb";
my $OUT    = "$BASE/out";

my $IMG_W = 3200;               # 3D image size in pixels. DS now draws NO text, so a big size is safe
my $IMG_H = 2400;
my $ZOOM  = 1.0;                # zoom after FitView so the molecule fills the picture
my $GFX_QUALITY = 'medium';       # low | medium | high | ultrahigh.  'ultrahigh' (the publication default) double-prints the text labels
my $TEST_QUALITY = 0;           # 1 = for the FIRST pdb also save one test picture per quality level (name_3D_q_<level>.png) -> choose the clean one
my $TEST_SIZES = 0;             # 1 = for the FIRST pdb also save 3 test sizes (to check label text); then set 0

my $LINE_WIDTH = 2;                # thickness of the residue lines (1 = very thin)
my $RECOLOUR_LABELS = 0;           # 1 = try to colour distance numbers (@DIST_COLOR); 0 = keep DS default (cleanest text)
my $LABEL_PUSH = 1.3;              # move each residue name this many Angstrom outward (away from the ligand) so the names spread out. 0 = off
my $SAVE_VIEWS = 1;                # 1 = also save extra viewing angles of every complex in out/views/ so you can pick the clearest
my @VIEWS = ([0,0,0], [90,0,0], [0,90,0], [45,45,0]);   # angles (degrees about x,y,z; each measured from the original orientation)
my @VIEW_ROT = (0, 0, 0);        # rotate the final view (degrees about x,y,z), e.g. (90,0,0) / (0,90,0) / (0,0,90) to choose the pose with least overlap
my @RESIDUE_COLOR = ();                 # () = element colours (O red, N blue, S yellow) with carbon in RESIDUE_CARBON; (r,g,b) = ONE colour for everything
my @DIST_COLOR = (200, 0, 0);           # colour of the distance numbers (dark red stands out from the dashed bond lines)
my @RESIDUE_CARBON = (0, 0, 0);        # carbon (and its bonds) of the interacting residues: black
my @LIGAND_CARBON = (0, 190, 190);    # ligand carbon colour (R,G,B). Cyan like your examples.
                                      # use (150,150,150) for grey carbon, (130,150,70) for olive
my $POCKET_EXTRA_RADIUS = 0;          # 0 = only interacting residues; e.g. 4.0 = also all residues within 4 A
my $LABEL_RESIDUES = 0;               # keep 0: names/distances are now written by finish_images.py (crisp text)
my $SHOW_DS_LINES = 0;                # 0 = DS hides its pale interaction lines; finish_images.py redraws them in strong colours. 1 = keep DS lines
my $CALIB_POINTS = 12;                # atoms used to learn the DS camera (2 marker pictures per view)
# ============================================================================

mkdir $BASE unless -d $BASE;
mkdir $OUT  unless -d $OUT;
open(my $LF, '>', "$OUT/batch_log.txt") or die "cannot write log in $OUT\n";
select((select($LF), $| = 1)[0]);
sub L { my $m = shift; print $LF "$m\n"; print "$m\n"; }

sub SetQuality {
    my ($doc, $q) = @_;
    if    ($q eq 'low')       { $doc->GraphicsQuality = Mdm::graphicsQualityLow; }
    elsif ($q eq 'medium')    { $doc->GraphicsQuality = Mdm::graphicsQualityMedium; }
    elsif ($q eq 'ultrahigh') { $doc->GraphicsQuality = Mdm::graphicsQualityUltrahigh; }
    else                      { $doc->GraphicsQuality = Mdm::graphicsQualityHigh; }
}

sub SiteInfo {
    my ($site) = @_;
    return ("", "") unless $site;
    my %rs;
    my $ats = eval { $site->Atoms };
    foreach my $at (@{ $ats || [] }) {
        my $p = eval { $at->Parent };
        my $rn = $p ? eval { $p->Name } : "";
        $rs{$rn} = 1 if $rn;
    }
    my $xyz = "";
    my $pt = eval { $site->XYZ };
    if ($pt) {
        my @c = map { my $v = eval { $pt->$_ }; defined $v ? sprintf("%.3f", $v) : "" } qw(X Y Z);
        $xyz = join(",", @c);
    }
    return (join(";", sort keys %rs), $xyz);
}

sub Try {                       # run code; log failure but never stop
    my ($name, $code) = @_;
    my @r = eval { $code->() };
    if ($@) { my $e = $@; $e =~ s/\s+$//; $e =~ s/\n.*//s; L("   FAILED $name :: $e"); return; }
    return (@r ? $r[0] : 1);
}

my @files = sort glob("$PDBDIR/*.pdb");
L("Found " . scalar(@files) . " pdb files in $PDBDIR");

my %METAL = map { $_ => 1 } qw(ZN MG CA FE MN CU CO NI NA K CD HG FE2 ZN2 MN2);

foreach my $pdb (@files) {
    my ($name) = $pdb =~ m{([^/]+)\.pdb$}i;
    L("");
    L("=== $name ===");
    my $ok = eval { ProcessOne($pdb, $name); 1 };
    L("   ABORTED $name :: $@") unless $ok;
}
L("");
L("ALL DONE");
close $LF;

# ------------------------------------------------------------------------------
sub ProcessOne {
    my ($pdb, $name) = @_;

    DiscoveryScript::Open($pdb);
    my $document = DiscoveryScript::LastActiveDocument(MdmModelType);
    die "cannot open $pdb\n" unless $document;
    $document->EnableUpdateViews(0);

    # ---- look: white background, quality, no cartoon ----
    SetPublicationQualityDisplay($document);
    Try("graphics quality $GFX_QUALITY", sub { SetQuality($document, $GFX_QUALITY); 1 });
    Try("hide protein ribbon", sub { $document->SetProteinDisplayStyle(Mdm::styleProteinNone); 1 });

    # ---- ligand / receptor ----
    my $ligands = GetLigands($document);
    my $prot    = GetProteinChains($document);
    die "no ligand found\n" unless $ligands && $ligands->Count > 0;
    L("   ligands: " . join(", ", map { $_->Name } @$ligands));

    # ---- interactions ----
    my $mon = eval { $document->CreateLigandNonbondMonitor({ RequirePredefinedLigand => 0 }) };
    die "CreateLigandNonbondMonitor failed: $@\n" unless $mon;
    Try("show all nonbonds + distances", sub {
        my $ok = eval { $mon->ShowAllNonbonds({ ShowDistanceLabels => 0 }); 1 };
        $mon->ShowAllNonbonds(1) unless $ok;
        1 });
    my $nbs = $document->Nonbonds;
    L("   interactions found: " . ($nbs ? $nbs->Count : 0));

    # ---- 2D diagram FIRST (needs the whole pocket still visible) ----
    my $f2 = "$OUT/${name}_2D.svg";
    Try("SaveTwoDView svg", sub { $document->SaveTwoDView($f2, 'svg'); 1 });
    L("   2D  -> " . (-s $f2 ? "$f2 (" . (-s $f2) . " bytes)" : "NOT CREATED"));

    # residues that REALLY interact (coloured circles in the 2D diagram; pale-green van der Waals circles are ignored)
    my %in2d;
    if (-s $f2 && open(my $sf, '<', $f2)) {
        local $/; my $svg = <$sf>; close $sf;
        my (@circ, @txt, %all);
        while ($svg =~ m{<circle\s+cx="([-\d.eE]+)"\s+cy="([-\d.eE]+)"[^>]*fill="url\(#(\w+)\)"}g) { push @circ, [$1,$2,$3]; }
        while ($svg =~ m{<text\s+x="([-\d.eE]+)"\s+y="([-\d.eE]+)"[^>]*>\s*([A-Za-z0-9]{1,4})\s*</text>\s*<text[^>]*>\s*([A-Za-z0-9]*):(-?\d+[A-Za-z]?)\s*</text>}g) {
            push @txt, [$1, $2, uc($3) . $5]; $all{uc($3) . $5} = 1;
        }
        foreach my $c (@circ) {
            next if $c->[2] eq 'rad_grad0';                     # van der Waals only -> skip
            my ($best, $bd) = (undef, 1e9);
            foreach my $t (@txt) {
                my $d = abs($t->[0] - $c->[0]) + abs($t->[1] - $c->[1]);
                ($best, $bd) = ($t, $d) if $d < $bd;
            }
            $in2d{$best->[2]} = 1 if $best && $bd < 1.0;
        }
        if (!%in2d) { %in2d = %all; L("   (could not separate vdW residues - using all contacts)"); }
    }
    L("   residues named in 2D diagram: " . scalar(keys %in2d));

    # ---- collect ligand + interacting residues (+ metals) ----
    my %ligId = map { $_->UniqueId => 1 } @$ligands;
    my (%seen, @show);
    foreach my $lg (@$ligands) { push @show, $lg; $seen{$lg->UniqueId} = 1; }

    foreach my $nb (@{ $nbs || [] }) {
        my $res = eval { $nb->Residues };
        foreach my $r (@{ $res || [] }) {
            next if $seen{$r->UniqueId}++;
            push @show, $r;
        }
    }
    my $all = eval { $document->Residues } || [];
    L("   residues in document: " . scalar(@$all));
    foreach my $r (@$all) {
        next if $seen{$r->UniqueId};
        my $nm = uc(eval { $r->Name } || '');
        if ($in2d{$nm}) { $seen{$r->UniqueId} = 1; push @show, $r; }
    }
    # metals close to ligand (Zn etc.)
    Try("metals near ligand", sub {
        $document->DeselectAll();
        foreach my $lg (@$ligands) { $lg->Select(); }
        $document->SelectByRadius(5.0, 'Atom');
        my $sel = $document->SelectedObjects();
        foreach my $a (@$sel) {
            next unless $a->IsA('Mdm::Atom');
            my $sym = uc(eval { $a->ElementSymbol } || '');
            next unless $METAL{$sym};
            my $p = $a->Parent;
            next if $seen{$p->UniqueId}++;
            push @show, $p;
        }
        $document->DeselectAll(); 1 });
    if ($POCKET_EXTRA_RADIUS > 0) {
        Try("extra pocket residues", sub {
            $document->DeselectAll();
            foreach my $lg (@$ligands) { $lg->Select(); }
            $document->SelectByRadius($POCKET_EXTRA_RADIUS, 'Atom');
            my $sel = $document->SelectedObjects();
            foreach my $a (@$sel) {
                next unless $a->IsA('Mdm::Atom');
                my $p = $a->Parent;
                next if $seen{$p->UniqueId}++;
                push @show, $p;
            }
            $document->DeselectAll(); 1 });
    }
    L("   residues kept in view: " . scalar(@show));

    # ---- table of interactions (type, distance, both ends) -> used by finish_images.py for the 2D distances ----
    Try("write nonbond table", sub {
        my $tsv = "$OUT/${name}_nonbonds.tsv";
        open(my $TF, '>', $tsv) or die "cannot write $tsv";
        print $TF "type\tdistance\tfrom_res\tfrom_xyz\tto_res\tto_xyz\n";
        my $n = 0;
        foreach my $nb (@{ $nbs || [] }) {
            my $d = eval { $nb->Distance };
            my $t = eval { $nb->TypeNames() };
            $t = "" unless defined $t;
            $t = join(";", @$t) if ref($t) eq 'ARRAY';
            my ($fr, $fx) = SiteInfo(eval { $nb->FromSite });
            my ($tr, $tx) = SiteInfo(eval { $nb->ToSite });
            printf $TF "%s\t%s\t%s\t%s\t%s\t%s\n", "$t", (defined $d ? sprintf("%.2f", $d) : ""), $fr, $fx, $tr, $tx;
            L("   nb: $t  d=" . (defined $d ? $d : "?") . "  from=[$fr] $fx  to=[$tr] $tx") if $n < 3;
            $n++;
        }
        close $TF;
        L("   nonbond table: $n rows -> $tsv");
        1 });

    # ---- what finish_images.py needs: ligand centre + one anchor point per residue ----
    my ($lcx, $lcy, $lcz, $nla) = (0, 0, 0, 0);
    foreach my $lg (@$ligands) {
        foreach my $at (@{ $lg->Atoms }) {
            my $p = eval { $at->XYZ } or next;
            $lcx += $p->X; $lcy += $p->Y; $lcz += $p->Z; $nla++;
        }
    }
    if ($nla) { $lcx /= $nla; $lcy /= $nla; $lcz /= $nla; }
    Try("write residue table", sub {
        open(my $RF, '>', "$OUT/${name}_labels.tsv") or die "cannot write labels table";
        print $RF "kind\ttext\tx\ty\tz\tcx\tcy\tcz\n";
        printf $RF "LIGAND\t\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $lcx, $lcy, $lcz, $lcx, $lcy, $lcz;
        my $n = 0;
        foreach my $o (@show) {
            next if $ligId{$o->UniqueId};
            my ($bx, $by, $bz, $bd) = (0, 0, 0, -1);
            my ($sx, $sy, $sz, $sn) = (0, 0, 0, 0);
            foreach my $at (@{ $o->Atoms }) {
                my $el = uc(eval { $at->ElementSymbol } || '');
                next if $el eq 'H';
                my $p = eval { $at->XYZ } or next;
                my ($x, $y, $z) = ($p->X, $p->Y, $p->Z);
                $sx += $x; $sy += $y; $sz += $z; $sn++;
                my $d = ($x-$lcx)**2 + ($y-$lcy)**2 + ($z-$lcz)**2;
                ($bx, $by, $bz, $bd) = ($x, $y, $z, $d) if $d > $bd;
            }
            next unless $sn;
            printf $RF "RES\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", $o->Name, $bx, $by, $bz, $sx/$sn, $sy/$sn, $sz/$sn;
            $n++;
        }
        close $RF;
        L("   residue table: $n residues -> ${name}_labels.tsv"); 1 });

    # ---- calibration atoms: well spread heavy atoms (farthest-point sampling) ----
    my @cand;
    foreach my $o (@show) {
        foreach my $at (@{ $o->Atoms }) {
            my $el = uc(eval { $at->ElementSymbol } || '');
            next if $el eq 'H';
            my $p = eval { $at->XYZ } or next;
            push @cand, [$at, $p->X, $p->Y, $p->Z];
        }
    }
    my @pts;     # [atom, pass, hue, x, y, z]
    if (@cand >= 8) {
        my ($mx, $my, $mz) = (0, 0, 0);
        foreach my $c (@cand) { $mx += $c->[1]; $my += $c->[2]; $mz += $c->[3]; }
        $mx /= @cand; $my /= @cand; $mz /= @cand;
        my ($first, $fd) = (undef, -1);
        foreach my $c (@cand) { my $d = ($c->[1]-$mx)**2 + ($c->[2]-$my)**2 + ($c->[3]-$mz)**2; ($first, $fd) = ($c, $d) if $d > $fd; }
        my @chosen = ($first);
        my %mind = map { ("$_" => 1e30) } @cand;
        while (@chosen < $CALIB_POINTS && @chosen < @cand) {
            my $last = $chosen[-1];
            my ($best, $bd) = (undef, -1);
            foreach my $c (@cand) {
                my $d = ($c->[1]-$last->[1])**2 + ($c->[2]-$last->[2])**2 + ($c->[3]-$last->[3])**2;
                $mind{"$c"} = $d if $d < $mind{"$c"};
                ($best, $bd) = ($c, $mind{"$c"}) if $mind{"$c"} > $bd;
            }
            last unless $best;
            push @chosen, $best;
        }
        for (my $i = 0; $i < @chosen; $i++) {
            my $c = $chosen[$i];
            push @pts, [$c->[0], $i % 2, int($i / 2) % 6, $c->[1], $c->[2], $c->[3]];
        }
    }
    Try("write calibration table", sub {
        open(my $CF, '>', "$OUT/${name}_calib.tsv") or die "cannot write calib table";
        print $CF "pass\thue\tx\ty\tz\n";
        printf $CF "%d\t%d\t%.3f\t%.3f\t%.3f\n", $_->[1], $_->[2], $_->[3], $_->[4], $_->[5] foreach @pts;
        close $CF;
        L("   calibration atoms: " . scalar(@pts)); 1 });

    my $ctx = { doc => $document, show => \@show, ligands => $ligands, ligId => \%ligId, mon => $mon, nbs => $nbs, pts => \@pts };

    # ---- look + view ----
    $document->EnableUpdateViews(1);
    ApplyLook($ctx);
    Try("fit view to pocket", sub {
        $document->DeselectAll();
        foreach my $o (@show) { $o->Select(); }
        $document->FitView();
        $document->DeselectAll(); 1 });
    Try("rotate view", sub { $document->RotateView(@VIEW_ROT, 0) if grep { $_ } @VIEW_ROT; 1 });
    Try("zoom", sub { $document->ScaleView($ZOOM); 1 }) if $ZOOM && $ZOOM != 1;
    $document->UpdateViews();

    mkdir "$OUT/calib" unless -d "$OUT/calib";
    Try("main picture", sub { Capture($ctx, $name, "v0", "$OUT/${name}_3D_geo.png"); 1 });
    L("   3D geometry -> " . (-s "$OUT/${name}_3D_geo.png" ? "OK" : "NOT CREATED"));

    if ($SAVE_VIEWS) {
        my $vdir = "$OUT/views";
        mkdir $vdir unless -d $vdir;
        my $vi = 0;
        foreach my $v (@VIEWS) {
            $vi++;
            Try("view $vi", sub {
                $document->RotateView($v->[0], $v->[1], $v->[2], 1);
                $document->UpdateViews();
                Capture($ctx, $name, "v$vi", "$vdir/${name}_view${vi}_geo.png");
                1 });
        }
        L("   extra views done (" . scalar(@VIEWS) . ")");
    }

    $document->Close();
}

# ------------------------------------------------------------------------------
# geometry picture + the marker pictures that tell finish_images.py where every 3D point lands on the screen
sub Capture {
    my ($ctx, $name, $tag, $geo) = @_;
    my $doc = $ctx->{doc};
    ApplyLook($ctx);
    $doc->UpdateViews();
    $doc->SaveImage($geo, 'png', $IMG_W, $IMG_H, 0);
    my @HUE = ([255,0,0], [255,255,0], [0,255,0], [0,255,255], [0,0,255], [255,0,255]);
    foreach my $pass (0, 1) {
        my @mine = grep { $_->[1] == $pass } @{ $ctx->{pts} };
        next unless @mine;
        my $ok = eval {
            $doc->DeselectAll();
            foreach my $p (@mine) { $p->[0]->Select(); }
            $doc->ShowOnly();
            foreach my $p (@mine) { eval { $p->[0]->SetVisible(); }; }
            $doc->SetAtomDisplayStyle({ DisplayStyle => 'styleAtomBallAndStick', StickSize => 0.3, BallSize => 0.9 });
            $doc->DeselectAll();
            foreach my $p (@mine) {
                $p->[0]->Select();
                $doc->SetColor(Mdm::Color::Create(@{ $HUE[$p->[2]] }));
                $doc->DeselectAll();
            }
            $doc->UpdateViews();
            $doc->SaveImage("$OUT/calib/${name}_${tag}_m$pass.png", 'png', $IMG_W, $IMG_H, 0);
            1 };
        L("   FAILED marker pass $pass ($tag) :: " . ($@ =~ /^(.*)/ ? $1 : '')) unless $ok;
    }
    ApplyLook($ctx);       # restore the normal picture
    $doc->UpdateViews();
}

# ------------------------------------------------------------------------------
# show only ligand + interacting residues, and give them the final style
sub ApplyLook {
    my ($ctx) = @_;
    my $document = $ctx->{doc};
    my @show = @{ $ctx->{show} };
    my $ligands = $ctx->{ligands};
    my %ligId = %{ $ctx->{ligId} };
    my $mon = $ctx->{mon};
    my $nbs = $ctx->{nbs};

    # ---- show only these ----
    $document->DeselectAll();
    foreach my $o (@show) { $o->Select(); }
    $document->ShowOnly();
    foreach my $o (@show) { $o->SetVisible(); }
    $document->DeselectAll();
    if ($SHOW_DS_LINES) {
        foreach my $nb (@{ $nbs || [] }) { eval { $nb->SetVisible(); }; }
        eval { $mon->SetVisible(); };
    }

    # ribbons of the visible residues come back when residues are shown -> hide them again
    my $chainList = eval { $document->AminoAcidChains };
    $chainList = [] unless $chainList;
    foreach my $ch (@$chainList) {
        eval { $ch->DisplayStyle = Mdm::styleProteinNone; };
    }
    eval { $document->SetProteinDisplayStyle(Mdm::styleProteinNone); };

    # ---- residues: thin lines, element colours, black carbons ----
    eval {
        $document->DeselectAll();
        foreach my $o (@show) { next if $ligId{$o->UniqueId}; $o->Select(); }
        $document->SetAtomDisplayStyle({ DisplayStyle => 'styleAtomLine' });
        $document->SetAtomColorScheme(Mdm::atomColorByElement);
        $document->DeselectAll(); };
    eval { $document->LineWidth = $LINE_WIDTH; };
    eval {
        $document->DeselectAll();
        foreach my $o (@show) {
            next if $ligId{$o->UniqueId};
            if (@RESIDUE_COLOR) { $o->Select(); }
            else { foreach my $a (@{ $o->Atoms }) { $a->Select() if (eval { $a->ElementSymbol } || '') eq 'C'; } }
        }
        $document->SetColor(Mdm::Color::Create(@RESIDUE_COLOR ? @RESIDUE_COLOR : @RESIDUE_CARBON));
        $document->DeselectAll(); };

    # ---- ligand: ball & stick, element colours, custom carbon colour ----
    eval {
        $document->DeselectAll();
        foreach my $lg (@$ligands) { $lg->Select(); }
        $document->SetAtomDisplayStyle({ DisplayStyle => 'styleAtomBallAndStick', StickSize => 0.28, BallSize => 0.32 });
        $document->SetAtomColorScheme(Mdm::atomColorByElement);
        $document->DeselectAll(); };
    eval {
        $document->DeselectAll();
        foreach my $lg (@$ligands) {
            foreach my $a (@{ $lg->Atoms }) { $a->Select() if (eval { $a->ElementSymbol } || '') eq 'C'; }
        }
        $document->SetColor(Mdm::Color::Create(@LIGAND_CARBON));
        $document->DeselectAll(); };
}
