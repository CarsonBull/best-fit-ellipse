#!/usr/bin/env perl

use strict;
use warnings;
use utf8;
use Test::More tests => 7;

my @TARGETS = ("efit", "eparam");
my @SRC = ("Makefile", "efit.c", "eparam.c");

foreach my $SRC (@SRC) {
	ok(-f "$SRC", "$SRC found");
}

my $CFLAGS="-Wall -pedantic";

foreach my $TARGET (@TARGETS) {
	my $compiler_output = `env CFLAGS="$CFLAGS" make -B $TARGET >/dev/null 2>/dev/null`;
	ok((!$? and -f "$TARGET"), "$TARGET built");
	$compiler_output = `env CFLAGS="$CFLAGS -Werror" make -B $TARGET 2>&1`;
	ok((!$? and -f "$TARGET"), "$TARGET built without warnings");
	if ($?) {
		$compiler_output =~ s/^/  /mg;
		diag "Compiler output:\n", $compiler_output;
	}
}
