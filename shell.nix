# This is a nix-shell for use with the nix package manager.
# If you have nix installed, you may simply run `nix-shell`
# in this repo, and have all dependencies ready in the new shell.

{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
  buildInputs = with pkgs;
    [
      gnumake
      clang-tools
      clang
      linuxKernel.packages.linux_6_6.perf  # performance tests

      # For experiments
      gnuplot
      timg
    ];
}
