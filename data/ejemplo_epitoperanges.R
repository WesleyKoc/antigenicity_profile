#Sa, Ca1, Ca2, Cb: Meng, F., Cheng, Z., Feng, Z., Zhang, Y., Zhang, Y., Wang, Y., Zhai, Y., Kuang, P., Qu, R., Chen, Y., Qiao, C., Chen, H., & Yang, H. (2025). H128N substitution in the SA antigenic site of HA1 causes antigenic drift between eurasian avian-like H1N1 and 2009 pandemic H1N1 influenza viruses. Viruses, 17(10), 1360. https://doi.org/10.3390/v17101360
#Sb: Matsuzaki, Y., Sugawara, K., Nakauchi, M., Takahashi, Y., Onodera, T., Tsunetsugu-Yokota, Y., Matsumura, T., Ato, M., Kobayashi, K., Shimotai, Y., Mizuta, K., Hongo, S., Tashiro, M., & Nobusawa, E. (2014). Epitope mapping of the hemagglutinin molecule of A/(H1N1)PDM09 influenza virus by using monoclonal antibody escape mutants. Journal of Virology, 88(21), 12364–12373. https://doi.org/10.1128/jvi.01381-14
H1N1_idx <- list("A"=c(128:129,158:167),#Technically Sa; Or 141:143,171:180
                 "B"=c(170:173,202:212),#Technically Sb; Or 170:173,202:212
                 "C"=c(169:173,206:208),#Technically Ca1; 
                 "D"=c(140:145,224:225),#Technically Ca2; Or 151:159,238
                 "E"=c(78:83))#Technically Cb

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