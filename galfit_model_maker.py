import os
import sys
import time
####
all_quasars = ['he1104','he2149','q2237']
####
print(' ')
print('Program to make galfit model needed for next galfitting!')
print('Currently, there are quasar models avaliable for:',all_quasars)
quasar = str(input('Input quasar name: '))
print(' ')
print('Your input was:',quasar)
####
if quasar not in all_quasars:
    print('Invalid quasar entry!')
    print(quasar,'is not in',all_quasars)
    sys.exit('Interrupting galfit_model_maker.py')
else:
    print('Valid entry!')
####
print('Starting to make',quasar,'in')
time.sleep(1)
print('3')
time.sleep(1)
print('2')
time.sleep(1)
print('1')
time.sleep(1)
####
print('Preparing',quasar+'.c, which takes approximately 42 seconds...')
os.system('cp saved_'+quasar+'_q2237.c q2237.c')
os.system('rm q2237.o')  
time.sleep(42)
####
print(' ')
print('Making "galfit" softlinked to all quasar folders!')
os.system('make')
time.sleep(1)
####
print(' ')
print('Check 5 last modified files in this folder:')
os.system('ls -ltr | tail -5')
print(' ')
time.sleep(1)
####
print('Galfit should now be ready for galfitting of',quasar+'!')
print(' ')
####
print('galfit_model_maker.py is finished!')
print(' ')