#Sa, Ca1, Ca2, Cb: Meng, F., Cheng, Z., Feng, Z., Zhang, Y., Zhang, Y., Wang, Y., Zhai, Y., Kuang, P., Qu, R., Chen, Y., Qiao, C., Chen, H., & Yang, H. (2025). H128N substitution in the SA antigenic site of HA1 causes antigenic drift between eurasian avian-like H1N1 and 2009 pandemic H1N1 influenza viruses. Viruses, 17(10), 1360. https://doi.org/10.3390/v17101360
#Sb: Matsuzaki, Y., Sugawara, K., Nakauchi, M., Takahashi, Y., Onodera, T., Tsunetsugu-Yokota, Y., Matsumura, T., Ato, M., Kobayashi, K., Shimotai, Y., Mizuta, K., Hongo, S., Tashiro, M., & Nobusawa, E. (2014). Epitope mapping of the hemagglutinin molecule of A/(H1N1)PDM09 influenza virus by using monoclonal antibody escape mutants. Journal of Virology, 88(21), 12364–12373. https://doi.org/10.1128/jvi.01381-14
H1N1_idx <- list("A"=c(141,142,143,162,163,171,172,173,174,177,180),#Technically Sa; Or 141:143,171:180
                 "B"=c(170,173,185,186,190,202,206,210,211,212),#Technically Sb; Or 170:173,202:212
                 "C"=c(186,187,188,189,190,220,222,223,224,225,237,252,253,255,256,257),#Technically Ca1; shifted by +17 in accordance to Meng et al from c(169, 170, 171, 172, 173,203,205,206,207,208,220,235,236,238,239,240)
                 "D"=c(137,138,141,151,154,156,157,158,159,200,222,238),#Technically Ca2; Or 151:159,238
                 "E"=c(90:94))#Technically Cb; shifted from 73:77 in accordance to Meng et al by +17

#Gupta V, Earl DJ, Deem MW. Quantifying influenza vaccine efficacy and antigenic distance. Vaccine. 2006 May 1;24(18):3881-8.
H3N2_idx <- list("A"=c(122,128,44,96,57,124,129,45,102,59,126,155,46,103,62,130,156,47,117,63,131,157,48,121,67),
                 "B"=c(132,158,50,167,75,133,159,51,170,78,135,160,53,171,80,137,163,54,172,81,138,165,273,173,82),
                 "C"=c(140,186,275,174,83,142,187,276,175,86,143,188,278,176,87,144,189,279,177,88,145,190,280,179,91),
                 "D"=c(146,192,294,182,92,150,193,297,201,94,152,194,299,203,109,168,196,300,207,260,197,304,208,261),
                 "E"=c(198,305,209,262,307,212,265,308,213,309,214,310,215,311,216,312,217,218,219,226,227,228,229,230,238,240,242,244,246,247,248))

#Pan, Y., & Deem, M. W. (2016). Prediction of influenza B vaccine effectiveness from sequence data. Vaccine, 34(38), 4610–4617. https://doi.org/10.1016/j.vaccine.2016.07.015
VIC_idx <- list("A"=c(121,122,123,125,126,134,135,136,137,139,141,142,144,146,147,148,149,150,151,155,157,177),
                "B"=c(127,129,133,160,161,162,163,164,165,166,168,172,174,196,197,198,199,200,202,203,204,206,207,208,209),
                "C"=c(34,35,36,37,38,39,40,289,291,292,293,294,309,315,317,318,320,321,323,324,325,326,327),
                "D"=c(93,101,102,116,120,176,179,180,182,183,184,185,186,187,188,190,212,214,218,219,220,223,224,225,226,
                    227,228,229,230,233,242,243,244,245,246,254,255,256,257,258),
                "E"=c(42,44,48,56,58,59,63,71,73,75,77,78,79,80,83,84,85,88,89,91,108,273,276,277,280))
