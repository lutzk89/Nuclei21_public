#include<iostream>
#include<list>
#include<fstream>
#include<cmath>
using namespace std;
class box {
    public:
        list<int> members;
        box new_box(list<int> members)
        {
            box* neu=new box;
            neu->members=members;
            return *neu;
        }
        int add_member(int member)
        {
            this->members.push_back(member);
            return 0;
        }
        int get_member()
        {
      	int output;
            if (members.size() > 0)
            {
                output = members.front();
                members.pop_front();
                return output;
            }
            else
                return -1;
        }
        list<int> get_members()
        {
            return this->members;
        }
        void clear()
        {
            this->members.clear();
        }
        int get_size()
        {
            return members.size();
        }
};
