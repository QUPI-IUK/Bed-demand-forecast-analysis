# startJuliaServer() before plan(multisession) ==>
# Environment variable will be picked up automatically
library(JuliaConnectoR)
library(future)
port <- startJuliaServer()
plan(multisession)
f1 <- future({
   #library(JuliaConnectoR)
   juliaEval("begin sleep(5); 1; end") # uses the Julia server running at the port
})

f2 <- future({
   #library(JuliaConnectoR)
   juliaEval("begin sleep(5); 2; end") # uses the Julia server running at the port
})

f1

f2

value(f1) # can be verified by the output here
value(f2)


#===============================================================================

library(JuliaConnectoR)
library(future)
plan(multisession)
port <- startJuliaServer()
f1 <- future({
   # to re-use the same Julia server, the environment variable must be set explicitly
   # if plan(multisession) is called before startJuliaServer().
   # Otherwise separate Julia servers are started by the child processes.
   Sys.setenv("JULIACONNECTOR_SERVER" = paste0("localhost:", port))
   juliaEval("begin sleep(5); 1; end")
})

f2 <- future({
   Sys.setenv("JULIACONNECTOR_SERVER" = paste0("localhost:", port))
   juliaEval("begin sleep(5); 2; end")
})

f1

f2

value(f1)
value(f2)


#===============================================================================
# Separate Julia processes for each R process

library(JuliaConnectoR)
library(future)
plan(multisession)
f1 <- future({
   juliaEval("begin sleep(5); 1; end")
})

f2 <- future({
   juliaEval("begin sleep(5); 2; end")
})

f1

f2

value(f1)
value(f2)



