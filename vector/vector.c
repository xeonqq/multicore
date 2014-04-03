#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;
int main(){
	vector<vector<int> > my_vector;
	int i,j;
	for(int i=0;i<10;i++){
		vector<int> sub_vector;
		for(int j=0;j<10;j++){
			sub_vector.push_back(i*10+j);
		}
		my_vector.push_back(sub_vector);
	}
	cout<<"the vector size is: " << my_vector.size()<<endl;
	vector<vector<int> >::iterator it;
	vector<int >::iterator it2;
	for(it = my_vector.begin(),i=0; it != my_vector.end(); ++it,i++)
	{
		for(it2 = (*it).begin(),j=0; it2 != (*it).end(); ++it2,j++)
		
		cout<<"the vector["<<i*my_vector[0].size()+j<<"] is "<<*it2 <<endl;
	}

	for( vector<int> int_vec : my_vector )
		for( int s : int_vec )
			cout << s << endl;
	
	return 0;
}
