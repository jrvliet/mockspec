       character*256 fname_gas, fname_part1, fname_part2, 
     &               fname_part3, fname_part4, 
     &               filename, sufix_name, middle_name, prefix_name 

       common /filenamesR/  fname_gas, fname_part1, fname_part2, 
     &               fname_part3, fname_part4, prefix_name

	 common /flags0/ iflag0a, iflag0b, iflag0c, 
     &           iflag0d, iflag0e, iflag0f
	 common /flags/ iflag1, iflag2, iflag3, iflag4, iflag5, 
     &           iflag6, iflag7, iflag8

	 common /GlobalV/ rho0C, gamma1, 
     &          rminC0, rmaxC0, ioptCenter, isys, Rs0,
     &          rmin0, rmax0, Zdmax0, nrbin, icomponent, 
     &          T_gas, tAge, lprefix
     
	 common /PreAna/ rminC,rmaxC, 
     &         xuser, yuser, zuser,Rs, AxX(3), AxY(3), AxZ(3)
	 common /Ana/ rmin,rmax,Zdmax



