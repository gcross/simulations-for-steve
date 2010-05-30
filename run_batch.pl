# foreach $number_of_sites (3,5,7,11,21,41,61,81,101) {
# foreach $lambda (0,0.1,0.2,0.4,0.8) {
foreach $number_of_sites (3,5,7) {
for ($lambda = 0; $lambda <= 4; $lambda += 0.1) {
    print "=====================================================================\n";
    system("./simulate",$lambda,$number_of_sites);
  }
}
