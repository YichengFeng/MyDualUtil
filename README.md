# Dual number for uncertainty propagation

Dual number is a way to calculate the precise number of (first-order) derivative of common math formula. 

Under normal C++ environment (e.g., c++ 9.4.0), the dual number classes are 
```bash
MyDualNumber.h # single variable
MyDualMultiv.h # multiple variables
```

Under cern ROOT environment (e.g., ROOT 6.18/05), I/O is realized by `IO\_ACLiC`, and visualization is realized by `TGraphErrors.h`, `TGraphAsymmErrors.h`.
```bash
MyDualGraph.h # statistical uncertainty
MySystGraph.h # statistical & systematic uncertainties
```

Each of those classes have their corresponding operators reloaded, which is the main idea of dual number.
