### SwatPestTileGw

Build with cmake:


### Configure the project

`FC=f77 FFLAGS="-ffree-line-length-0 -ffixed-line-length-0 -g -O2" cmake -B build`

### Build

`cmake --build build`

The binary `swat_681` will be copied to the root directory of this project folder. 

### Rebuild
Delete the `build` directory and `swat_681` binary then run the above commands again. Or run:

`cmake --build build --clean-first`
