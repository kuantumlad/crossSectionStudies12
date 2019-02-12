{
  ifstream aaa("clas12.lund");
  float px, py, pz, ex, ey, ez;
  float dummy;
  float phiE, phiP;
  for (int l = 0; l < 20; l++){
    for (int t = 0; t < 38; t++){
      aaa >> dummy;
      if (t == 16) ex = dummy;
      if (t == 17) ey = dummy;
      if (t == 18) ez = dummy;

      if (t == 30) px = dummy;
      if (t == 31) py = dummy;
      if (t == 32) pz = dummy;
    }
    phiE = atan2(ey, ex);
    phiP = atan2(py, px);
    cout << "phiE: " << phiE*57 << endl;
    cout << "phiP: " << phiP*57 << endl;
    cout << "px: " << px << endl;
    cout << "py: " << py << endl;
    cout << "ex: " << ex << endl;
    cout << "ey: " << ey << endl;

  }
}
