#define elDefaultPIDCut_hh

bool ele_default_PID_cut(int j){
  if(vpart_pid->at(j) == 11) return true;
  else return false;
}
