# var mysession1;
# var myage =30;
# ocpu.call("listQx",{age:myage,numberSims :100} , 
#           // first call 
#           function(session1){
#             mysession1 = session1;
#             session1.getConsole(function(txt){alert(session1.getKey())});
#             myKey = session1.getKey();
#             alert(JSON.stringify(myKey))
#             alert(myKey)
#             ocpu.call("drawLifeExpectancy",{qxStochasticSet:myKey,age:30},
#                       function(session2){
#                         mysession2 = session2;
#                         alert(session2.getKey())
#                       }
#             );
#             
#           }
#           
#           
# )
# 
# 
# alert(mysession1.getKey())