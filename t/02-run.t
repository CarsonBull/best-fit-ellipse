#!/usr/bin/env perl

use strict;
use warnings;
use utf8;

use Test::More tests => 4;

my $number_regex = qr/(?:[+-])?(nan|inf|\d+(?:\.\d+)(?:e[+-]\d+)?)/;

ok(-f "efit", "efit binary exists") || die "No point continuing with no binary";

my $efit_input = "5.9179400456751 1.06414138471654\n5.85063278484678 1.63604670195337\n5.7025094282751 2.31525295163843\n5.35903328637879 3.01059071850842\n4.82436358685989 3.58547989134214\n4.08591673692726 4.10129241378114\n3.3278601145453 4.49869655416849\n2.62702720808155 4.8650836520863\n1.72845159056454 4.90266694486781\n0.864359097691146 4.91535129136555\n-0.184037825609909 4.94438766901306\n-0.875838886311746 4.74519441278574\n-1.68796397397809 4.26662897873927\n-2.40980577607129 3.92480985418556\n-3.13761880365293 3.34518992006853\n-3.4489255745481 2.82273236273401\n-3.72886014027258 2.14752805258506\n-4.0378946975367 1.45749347547436\n-4.06242531269368 0.55026706057986\n-3.7451598211439 0.0243768134542408\n-3.55533912109331 -0.76618443688855\n-2.98131041800736 -1.36515442854396\n-2.39881663423544 -1.79171974232335\n-1.6641376566391 -2.24443911723967\n-1.05470301541246 -2.69948280223873\n-0.20788600304893 -2.98047829691554\n0.785603159047683 -3.06480407686772\n1.77013125061392 -2.98826639318711\n2.55864936066606 -2.78601525863493\n3.30948285276936 -2.50787586781035\n4.20799588546615 -2.11959444141323\n4.67709067551136 -1.57649682413335\n5.19383639936424 -1.04927888005981\n5.71587726734955 -0.38528615285236\n5.99339712780838 0.193725735264042\n";
my $efit_output = `printf '%s' "$efit_input" | ./efit`;
if ($?) {
	fail("efit crashed");
}
else {
	my @efit_nums = ($efit_output =~ /$number_regex/g);
	my $efit_error = pop @efit_nums;

	cmp_ok(abs($efit_error), '<', 1, "efit error < 1");
}

ok(-f "eparam", "eparam binary exists") || die "No point continuing with no binary";
my $eparam_input = "0.0056119275\n1.5463826460\n-1.9877179482\n-3.0868914114\n-22.2769694769\niters=257,error=0.742130396007\n";
my $eparam_output = `printf '%s' "$eparam_input" | ./eparam`;
if ($?) {
	fail("eparam crashed");
}
else {
	my @eparam_nums = ($eparam_output =~ /$number_regex/g);
	my $eparam_error = pop @eparam_nums;

	cmp_ok(abs($eparam_error), '<', 1, "eparam error < 1");
}
