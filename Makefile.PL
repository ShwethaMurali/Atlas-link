use 5.008008;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

WriteMakefile(
    NAME              => 'Atlas::vEXAMer',
    #VERSION_FROM      => 'src/From_Atlas/vEXAMer.pm', # finds $VERSION
    VERSION_FROM      => 'lib/Atlas/vEXAMer.pm', # finds $VERSION
    PREREQ_PM         => {
						 Statistics::Descriptive::Weighted => 0,
						 Graph => 0,
						 Statistics::Descriptive => 0,
						 XML::DOM => 0,
						 XML::Parser => 0,
						 XML::RegExp => 0,
						 Algorithm::ClusterPoints => 0
						 }, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Atlas/vEXAMer.pm', # retrieve abstract from module
       AUTHOR         => 'Jixin Deng <jdeng@bcm.edu>') : ()),
);
