all: PointPolygon.so u5m.so
    # Comment here preserves the prior tab
PointPolygon.so: PointPolygon.cpp
    Rscript --vanilla -e "TMB::compile('PointPolygon.cpp','-O0 -g')"
u5m.so: u5m.cpp
    Rscript --vanilla -e "TMB::compile('u5m.cpp','-O0 -g')"
