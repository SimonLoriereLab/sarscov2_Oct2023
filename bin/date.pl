use DateTime;

my $start_of_week =
    DateTime->new(
	year => 2023,
	month => 05,
	day => 25,
    )->truncate( to => 'week' );

print($start_of_week."\n");
