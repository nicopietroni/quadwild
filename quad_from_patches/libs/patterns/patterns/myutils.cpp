#include "myutils.h"

namespace utility{


    bool ispermuted(const vector<int>& list1,const vector<int>& list2, vector<pair<int,int>>& correspondence ){
        size_t n=list1.size();
        assert(list2.size()==n);
        correspondence.resize(n);
        int target=list1[0];
        for(size_t i=0;i<n;i++){
            if(list2[i]==target){
                bool yes=true;
                for(size_t j=0;j<n;j++){
                    correspondence[j]=make_pair(j,(i+j)%n);
                    yes=yes && (list2[(i+j)%n]==list1[j]);
                    if(yes==false)
                        break;
                }
                if(yes==true)
                   return true;
                yes=true;
                for(size_t j=0;j<n;j++){
                    correspondence[j]=make_pair(j,(i-j+n)%n);
                    yes=yes && (list2[(i-j+n)%n]==list1[j]);
                    if(yes==false)
                        break;
                }
                if(yes==true)
                   return true;
            }
        }
        return false;
    }

}
