#include <bits/stdc++.h>
using namespace std;

void separate_el2(int np, int ne, int nd[][2], int* ncrv, int parts[], int indices[]){
  // 添字を0からにする
  for(int i=0; i<ne; i++){
    nd[i][0]--;
    nd[i][1]--;
  }

  // a[i] = i番nodeを始点とするedgeの番号 O(n)
  vector<int> a(np);
  for(int i=0; i<ne; i++){
    a[nd[i][0]] = i;
  }

  // succ[i] = i番edgeの次のedgeの番号 O(n)
  vector<int> succ(ne);
  for(int i=0; i<ne; i++){
    succ[i] = a[nd[i][1]];
  }

  // edgeの番号の集まり
  // unordered_set<int> s;
  set<int> s;
  for(int i=0; i<ne; i++){
    s.insert(i);
  }  

  vector<vector<int>> curves;
  
    
  while(!s.empty()){
    vector<int> curve;
    
    // sからpop
    int start = *s.begin();
    s.erase(s.begin());
    // curveに追加
    curve.emplace_back(start);

    while(succ[*curve.rbegin()] != *curve.begin()){
      int v = succ[*curve.rbegin()];
      curve.emplace_back(v);
      
      s.erase(v);     				    
      
      if(s.empty()) break;
      
    }
    
    curves.emplace_back(curve);

  }

  *ncrv = curves.size();

  int i = 0;
  // for(auto it=curves.begin(); it!=curves.end(); it++){
  for(int icrv=0; icrv<curves.size(); icrv++){
    parts[icrv] = curves[icrv].size();
    for(int j=0; j<curves[icrv].size(); j++){
      indices[i] = curves[icrv][j] + 1; // 添字は1から
      i++;
    }
  }

  
  
}


void gen_curves(int np, int ne, int nd[][2], vector<vector<int>> *curves){
  // 添字を0からにする
  for(int i=0; i<ne; i++){
    nd[i][0]--;
    nd[i][1]--;
  }

  // a[i] = i番nodeを始点とするedgeの番号 O(n)
  vector<int> a(np);
  for(int i=0; i<ne; i++){
    a[nd[i][0]] = i;
  }

  // succ[i] = i番edgeの次のedgeの番号 O(n)
  vector<int> succ(ne);
  for(int i=0; i<ne; i++){
    succ[i] = a[nd[i][1]];
  }

  // edgeの番号の集まり
  // unordered_set<int> s;
  set<int> s;
  for(int i=0; i<ne; i++){
    s.insert(i);
  }  

  // vector<vector<int>> curves;
    
  while(!s.empty()){
    vector<int> curve;
    
    // sからpop
    int start = *s.begin();
    s.erase(s.begin());
    // curveに追加
    curve.emplace_back(start);

    while(succ[*curve.rbegin()] != *curve.begin()){
      int v = succ[*curve.rbegin()];
      curve.emplace_back(v);
      
      s.erase(v);     				    
      
      if(s.empty()) break;
      
    }
    
    (*curves).emplace_back(curve);

  }
}

int main(){
  // el2を標準入力から読み込み
  int np, tmp;
  scanf("%d", &np);
  double p[np][2];
  for(int i=0; i<np; i++){
    cin >> tmp >> p[i][0] >> p[i][1];   
  }
  int ne;
  scanf("%d", &ne);
  int nd[ne][2];
  // int nd;
  for(int i=0; i<ne; i++){
    cin >> tmp >> nd[i][0] >> nd[i][1];
  }

  int ncrv;
  int parts[ne];
  int indices[ne];
  separate_el2(np, ne, nd, &ncrv, parts, indices);

  cout << ncrv << endl;

  int j=0;
  for(int icrv=0; icrv<ncrv; icrv++){
    for(int i=0; i<parts[icrv]; i++){
      cout << indices[j] << " ";
      j++;
    }
    cout << endl;
  }

  // vector<vector<int>> curves;

  // gen_curves(np, ne, nd, &curves);

  // for(auto it=curves.begin(); it!=curves.end(); it++){
  //   cout << (*it).size() << endl;
  // }
  
}

