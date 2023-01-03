# This is a nix-shell for use with the nix package manager.
# If you have nix installed, you may simply run `nix-shell`
# in this repo, and have all dependencies ready in the new shell

{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
  buildInputs = with pkgs;
    [
      gnumake
      clang-tools
      clang
      fftw
      linuxKernel.packages.linux_5_15.perf  # performance tests
      cudaPackages_11_8.cudatoolkit 

      # For experiments
      gnuplot
      timg
    ];
  shellHook = ''
#    export CUDA_PATH=${pkgs.cudaPackages_11_8.cudatoolkit}
#    export LD_LIBRARY_PATH=${pkgs.cudaPackages_11_8.cudatoolkit}/lib:$LD_LIBRARY_PATH
#    export CPATH=${pkgs.cudaPackages_11_8.cudatoolkit}/include:$CPATH
   '';
}
