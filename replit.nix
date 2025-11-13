{ pkgs }: {
  deps = [
    pkgs.expat
    pkgs.python311
    pkgs.python311Packages.pip
    pkgs.python311Packages.flask
  ];
}
