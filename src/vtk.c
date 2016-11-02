
void WriteVTK(VecI box,int framenum){
	const int  ell = 0;
	const int disk = 1;
	double       dtype;
	VecR cmass; //wrapped Centre of mass
	// variables for colouring
	double X2,Y2,Z2;
	double YoX, ZoY; // Y2/X2 , Z2/Y2
	double RR;

	std::stringstream filename;
	filename<<"ellipsoids_N"<<numofDendrimers<<"_b"<<box.x<<".vtk";
	ofstream outf(filename.str().c_str(),std::ios::out);
	outf<<"# vtk DataFile Version 3.0"<<endl;
	outf<<"Random data to test tensors"<<endl;
	outf<<"ASCII"<<endl;
	outf<<"DATASET UNSTRUCTURED_GRID"<<endl;
	outf<<"POINTS "<<numofDendrimers<<" float"<<endl;
  	for (int dendid = 0;dendid<numofDendrimers;dendid++){
		cmass = D.at(dendid).getCMassFrame(framenum);			
		PBCAll(cmass,box);
		outf<<cmass.x<<" "<<cmass.y<<" "<<cmass.z<<endl;
	}
	outf<<endl;
	
	outf<<"POINT_DATA "<<numofDendrimers<<endl;
	outf<<"SCALARS scalars1 float"<<endl;
	outf<<"LOOKUP_TABLE default"<<endl;
  	for (int dendid = 0;dendid<numofDendrimers;dendid++) outf<<dendid<<" ";
	outf<<endl<<endl;
	outf<<"SCALARS RR float"<<endl;
	outf<<"LOOKUP_TABLE default"<<endl;
	//scale by eigenvalues
  	for (int dendid = 0;dendid<numofDendrimers;dendid++){
		X2   = D.at(dendid).getdblFrame("X2", framenum);
		Y2   = D.at(dendid).getdblFrame("Y2", framenum);
		Z2   = D.at(dendid).getdblFrame("Z2", framenum);
		YoX  = Y2 / X2;
		ZoY  = Z2 / Y2;

		RR   = ( Y2 - ( X2 + Z2)/2 )/( ( X2 - Z2)/2 );

		outf<<RR<<" ";
	}
	outf<<endl<<endl;
	outf<<"TENSORS tensors1 float"<<endl;
  	for (int dendid = 0;dendid<numofDendrimers;dendid++){
		std::vector<double> gyrtensor = D.at(dendid).getGyrTenFrame(framenum);
		std::vector<double>::iterator it;
		int i=0;
		for(it=(gyrtensor).begin();it!=(gyrtensor).end();++it){
			i++;
			outf<<*it<<" ";
			if (i==3) {outf<<endl;i=0;}
		}
		outf<<endl;
	}
	outf.close();
}
