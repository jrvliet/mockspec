import numpy as np
from pylab import *

## blzbinc,blzmean,blzstddn,blzstdup = np.genfromtxt('2blueloz_wrtzabs.tpcf',
##                                                   unpack=True)
## bhzbinc,bhzmean,bhzstddn,bhzstdup = np.genfromtxt('2bluehiz_wrtzabs.tpcf',
##                                                   unpack=True)
rlzbinc,rlzmean,rlzstddn,rlzstdup = np.genfromtxt('2redloz_wrtzabs.tpcf',
                                                  unpack=True)
rhzbinc,rhzmean,rhzstddn,rhzstdup = np.genfromtxt('2redhiz_wrtzabs.tpcf',
                                                  unpack=True)



#fig,((p11,p12),(p21,p22)) = subplots(2,2,figsize=(10.2,9))
fig,p22 = subplots(1,1,figsize=(5,4.5))

## p11.fill_between(blzbinc,blzstdup,blzstddn,color='#4169E1',alpha=0.3)
## p11.fill_between(rlzbinc,rlzstdup,rlzstddn,color='#DC143C',alpha=0.2)

## p11bl, = p11.plot(blzbinc,blzmean,color='#4169E1',lw=1.5)
## p11rd, = p11.plot(rlzbinc,rlzmean,color='#DC143C',lw=1.4)


## p12.fill_between(bhzbinc,bhzstdup,bhzstddn,color='#4169E1',alpha=0.3)
## p12.fill_between(rhzbinc,rhzstdup,rhzstddn,color='#DC143C',alpha=0.2)

## p12.plot(bhzbinc,bhzmean,color='#4169E1',lw=1.5)
## p12.plot(rhzbinc,rhzmean,color='#DC143C',lw=1.4)


## p21.fill_between(blzbinc,blzstdup,blzstddn,color='#32CD32',alpha=0.3)
## p21.fill_between(bhzbinc,bhzstdup,bhzstddn,color='#800080',alpha=0.2)

## p21lz, = p21.plot(blzbinc,blzmean,color='#228F22',lw=1.5)
## p21hz, = p21.plot(bhzbinc,bhzmean,color='#800080',lw=1.4)


p22.fill_between(rlzbinc,rlzstdup,rlzstddn,color='#32CD32',alpha=0.3)
p22.fill_between(rhzbinc,rhzstdup,rhzstddn,color='#800080',alpha=0.2)

p22lz, = p22.plot(rlzbinc,rlzmean,color='#228F22',lw=1.5)
p22hz, = p22.plot(rhzbinc,rhzmean,color='#800080',lw=1.4)



## p11leg = p11.legend([p11bl,p11rd],[r'Blue}',r'Red'],loc=0)
## p12leg = p12.legend([p11bl,p11rd],[r'Blue}',r'Red'],loc=0)
## p21leg = p21.legend([p21lz,p21hz],
##                     [r'Low $z_{gal}$',
##                      r'High $z_{gal}$'],loc=0)
p22leg = p22.legend([p22lz,p22hz],
                    [r'Low $z_{gal}$',
                     r'High $z_{gal}$'],loc=0)

for a in [p22leg]: #[p11leg,p12leg,p21leg,p22leg]:
    ltext = a.get_texts()
    setp(ltext,fontsize=11)

[a.set_xlim([0,350]) for a in [p22]] #[p11,p12,p21,p22]]
[a.set_ylim([0,0.14]) for a in [p22]] #[p11,p12,p21,p22]]

## p11.text(200,0.08,r'Low $z_{gal}$')
## p12.text(200,0.08,r'High $z_{gal}$')

## p21.text(200,0.08,r'Blue')
p22.text(200,0.08,r'Red')

[a.set_xlabel(r'$\Delta v_{cloud}$') for a in [p22]]
# [p11,p12,p21,p22]]

[a.set_ylabel(r'P$(\Delta v_{cloud})$') for a in [p22]]
 #[p11,p12,p21,p22]]

subplots_adjust(hspace=0.21,wspace=0.28)

savefig('2BKz_wrtzabs.pdf',bbox_inches='tight')
#savefig('2BKz_wrtzabs.eps',bbox_inches='tight')
