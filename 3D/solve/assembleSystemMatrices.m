function [H, G, HO, GO, HIn, GIn] = assembleSystemMatrices(HElaBou, GElaBou, ...
    HOrig, GOrig, HElaIn, GElaIn, tranU, tranT)
H = HElaBou/tranU;
G = GElaBou/tranT;
HO = HOrig/tranU;
GO = GOrig/tranT;
HIn = HElaIn/tranU;
GIn = GElaIn/tranT;
end