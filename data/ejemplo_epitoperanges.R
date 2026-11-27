#Deem MW, Pan K. The epitope regions of H1-subtype influenza A, with application to vaccine efficacy. Protein Eng Des Sel. 2009 Sep;22(9):543-6.
#https://mwdeem.org/pepitope/
H1N1_idx <- list("A"=c(118,120,121,122,126,127,128,129,132,133,134,135,137,139,140,141,142,143,146,147,149,165,252,253),
                 "B"=c(124,125,152,153,154,155,156,157,160,162,183,184,185,186,187,189,190,191,193,194,195,196),
                 "C"=c(34,35,36,37,38,40,41,43,44,45,269,270,271,272,273,274,276,277,278,283,288,292,295,297,298,302,303,305,306,307,308,309,310),
                 "D"=c(89,94,95,96,113,117,163,164,166,167,168,169,170,171,172,173,174,176,179,198,200,202,204,205,206,207,208,209,210,211,212,213,214,215,216,222,223,224,225,226,227,235,237,239,241,243,244,245),
                 "E"=c(47,48,50,51,53,54,56,57,58,66,68,69,70,71,72,73,74,75,78,79,80,82,83,84,85,86,102,257,258,259,260,261,263,267))

#Sa, Ca1, Ca2, Cb: Meng, F., Cheng, Z., Feng, Z., Zhang, Y., Zhang, Y., Wang, Y., Zhai, Y., Kuang, P., Qu, R., Chen, Y., Qiao, C., Chen, H., & Yang, H. (2025). H128N substitution in the SA antigenic site of HA1 causes antigenic drift between eurasian avian-like H1N1 and 2009 pandemic H1N1 influenza viruses. Viruses, 17(10), 1360. https://doi.org/10.3390/v17101360
#Sb: Matsuzaki, Y., Sugawara, K., Nakauchi, M., Takahashi, Y., Onodera, T., Tsunetsugu-Yokota, Y., Matsumura, T., Ato, M., Kobayashi, K., Shimotai, Y., Mizuta, K., Hongo, S., Tashiro, M., & Nobusawa, E. (2014). Epitope mapping of the hemagglutinin molecule of A/(H1N1)PDM09 influenza virus by using monoclonal antibody escape mutants. Journal of Virology, 88(21), 12364–12373. https://doi.org/10.1128/jvi.01381-14
H1N1_idx <- list("A"=c(141,142,143,162,163,171,172,173,174,177,180),#Technically Sa; Or 141:143,171:180
                 "B"=c(170,173,185,186,190,202,206,210,211,212),#Technically Sb; Or 170:173,202:212
                 "C"=c(186,187,188,189,190,220,222,223,224,225,237,252,253,255,256,257),#Technically Ca1; shifted by +17 in accordance to Meng et al from c(169, 170, 171, 172, 173,203,205,206,207,208,220,235,236,238,239,240)
                 "D"=c(137,138,141,151,154,156,157,158,159,200,222,238),#Technically Ca2; Or 151:159,238
                 "E"=c(90:94))#Technically Cb; shifted from 73:77 in accordance to Meng et al by +17

#A, B, D: Thompson, A. J., Wu, N. C., Canales, A., Kikuchi, C., Zhu, X., de Toro, B. F., Cañada, F. J., Worth, C., Wang, S., McBride, R., Peng, W., Nycholat, C. M., Jiménez-Barbero, J., Wilson, I. A., & Paulson, J. C. (2024). Evolution of human H3N2 influenza virus receptor specificity has substantially expanded the receptor-binding domain site. Cell Host & Microbe, 32(2), 261–275.e4. https://doi.org/10.1016/j.chom.2024.01.003
#C: Broecker, F., Sean, Sun, W., Krammer, F., Simon, V., & Palese, P. (2018). Immunodominance of antigenic site B in the hemagglutinin of the current H3N2 influenza virus in humans and mice. Journal of Virology, 92(20). https://doi.org/10.1128/jvi.01100-18
#E: Yamayoshi, S., Hamabata, T., Inoue, Y., Kiso, M., Okuda, M., Akasaka, H., Shihoya, W., Imamura, T., Matsuzawa, Y., Jubishi, D., Iwatsuki-Horimoto, K., Wilson, P. C., Nureki, O., & Kawaoka, Y. (2026). Characterisation of a human monoclonal antibody targeting a conserved epitope at the base of the HA head of influenza A(H3N2) virus. eBioMedicine, 131, Article 106450. https://doi.org/10.1016/j.ebiom.2026.106450
H3N2_idx <- list("A"=c(122:146),#critical cluster-transition positions 135, 144, 145
                 "B"=c(155:163,186:197),#immunodominant residues 156, 158, 159, 189, 193
                 "C"=c(45:54,275:278),
                 "D"=c(201:228),
                 "E"=c(62:83))

#A, C, D, E: Wan, X.-F., Guan, M., Balamalaliyage, P., Chen, H., Prasai, K., Alcala, A., Driver, J., Olivier, A. K., Gu, W., Frymire, C., Madrid, De D. M. C., Gao, C., Wang, C., Li, T., Tunterak, W., Yang, Q., Ramesh, A., Ali, M., Smith, D., Li, L., Sant, A. J., Hang, J., Xie, H., Zhou, M., Tao, Y. J. (2026). Epitope-spanning antigenic variation reprograms immunodominance and broadens immunity in sequential influenza vaccination. Nature Communications, 17. https://doi.org/10.1038/s41467-026-70202-y
#B: Guo, C.-Y., Zhang, H.-X., Zhang, J.-J., Sun, L.-J., Li, H.-J., Liang, D.-Y., Feng, Q., Li, Y., Feng, Y.-M., Xie, X., & Hu, J. (2019). Localization analysis of heterophilic antigen epitopes of H1N1 influenza virus hemagglutinin. Virologica Sinica, 34(3), 306–314. https://doi.org/10.1007/s12250-019-00100-9
VIC_idx <- list("A"=c(122,126,130:140,145:150),
                "B"=c(155:160,166,209,186:199),
                "C"=c(40,280:300),
                "D"=c(233,245:255),
                "E"=c(62:63,73,75:82,85:100))


if (current_subtype() == "H1N1") {
  A<-c(128:129,158:167)#Technically Sa; Or 141:143,171:180
  B<-c(170:173,202:212)#Technically Sb; Or 170:173,202:212
  C<-c(169:173,206:208)#Technically Ca1; 
  D<-c(140:145,224:225)#Technically Ca2; Or 151:159,238
  E<-c(78:83)#Technically Cb
} else if (current_subtype()  == "H3N2") {
  A<-c(122:146)#critical cluster-transition positions 135, 144, 145
  B<-c(155:163,186:197)#immunodominant residues 156, 158, 159, 189, 193
  C<-c(45:54,275:278)
  D<-c(201:228)
  E<-c(62:83)
} else {
  A<-c(122,126,130:140,145:150)
  B<-c(155:160,166,209,186:199)
  C<-c(40,280:300)
  D<-c(233,245:255)
  E<-c(62:63,73,75:82,85:100)
}




#New H1N1 epitopes 10.1016/j.vaccine.2012.07.079
if (current_subtype() == "H1N1") {
  A<-c(121, 127, 128, 130, 133, 141, 146)#Epitope Site A
  B<-c(125, 153, 160, 183, 186, 189, 190, 191, 193, 194)#Epitope Site B
  C<-c(35, 36, 43, 271, 273, 274, 277, 295, 310)#Epitope Site C
  D<-c(94, 163, 205, 209, 216, 222, 224)#Epitope Site D
  E<-c(47, 54, 69, 71, 73, 80, 82, 267)#Epitope Site E
}
#Deem MW, Pan K. The epitope regions of H1-subtype influenza A, with application to vaccine efficacy. Protein Eng Des Sel. 2009 Sep;22(9):543-6.
#https://mwdeem.org/pepitope/
H1N1_idx <- list("A"=c(121,127,128,130,133,141,146),
                 "B"=c(125,153,160,183,186,189,190,191,193,194),
                 "C"=c(35,36,43,271,273,274,277,295,310),
                 "D"=c(94,163,205,209,216,222,224),
                 "E"=c(47,54,69,71,73,80,82,267))