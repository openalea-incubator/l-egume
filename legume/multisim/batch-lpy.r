#tets lancement en R du batch l-py

dir <- choose.dir()
##"C:\\devel\\l-egume\\legume"

setwd(dir)



#single run avec options
system('python l-egume_run.py -u liste_usms_mix.xls SimTest 3')

#uilisable pour faire de l'optimisation? -> a tester


##"C:\\devel\\l-egume\\legume\\multisim\\sorties"
setwd(dir)
system('python l-egume_batch_mixture8.py')



#https://stackoverflow.com/questions/34449241/running-windows-command-prompt-using-r
#https://stackoverflow.com/questions/10155703/call-python-with-system-in-r-to-run-a-python-script-emulating-the-python-conso
