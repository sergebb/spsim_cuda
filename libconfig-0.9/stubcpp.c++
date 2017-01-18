
#include <stdio.h>
#include "libconfig.h++"

using namespace libconfig;

int main(int argc, char **argv)
  {
  Config cfg;

  try
    {
    FILE *fp = fopen("test.cfg", "r");
    cfg.read(fp);
    fclose(fp);
    
//  const ConfigSetting& setting = cfg.lookup("application.window.size.w");
//  long val = setting;
    
    long val = cfg.lookup("application.window.size.w");
    printf("val: %ld\n", val);
    
    Setting &s = cfg.lookup("application.group1.my_array");  
    long val4 = s[4];
    printf("item #4 is: %d\n", val4);

    Setting &grp = cfg.lookup("application.group1.group2");

    Setting &root = cfg.getRoot();

    Setting &rootn = root.add("new-one-at-top", TypeGroup);

    Setting &ngp = rootn.add("element", TypeFloat);
    
    ngp = 1.1234567890123;
    
//    long val22 = s[22];
//    printf("item #22 is: %d\n", val22);

    Setting &snew = grp.add("foobar", TypeArray);

    snew.add(TypeInt);
    snew.add(TypeInt);

    snew.add(TypeInt);
    snew.add(TypeInt);
    
    puts("created new array");

    snew[0] = 55;
    puts("elem 0");
    snew[1] = 66;
    puts("elem 1");

    Setting &sdel = cfg.lookup("application");

//    sdel.remove("group1");

    
    

    cfg.write(stdout);

    
    }
  catch(ParseException& ex)
    {
    printf("error on line %d: %s\n", ex.getLine(),
           ex.getError());
    }
    
  return(0);
  }
