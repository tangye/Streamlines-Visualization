#include"myTools.h"
	//
void writePara(char* strs,double para1,double para2 ,double para3,double para4)
{
	ofstream vtkfile;
	vtkfile.open ("Paraout.out", ios::out | ios::app);
	if(vtkfile.is_open())
	{
		vtkfile<<strs<<"\t"<<para1<<"\t"<<para2<<"\t"<<para3<<"\t"<<para4<<endl;
	}
	vtkfile.close();
};
	//