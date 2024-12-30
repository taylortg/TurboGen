The goal with the backend is to have a stand alone implementation of the model that is not gui driven
That way the model can get developed and abstracted so the gui only relies on the abstraction.
To-do:
[ ] Add in some terminal controls
[x] Finish implementing single zone Japikse model
[ ] Review radcomp
[ ] Implement the radcomp model
[ ] Review other github repo found (forgot the name)
** Both github repos rely on Aungier and Galvas models. The crux of the models are to calculate ideal
   impeller outlet conditions. Then use efficiency or enthalpy decrements to get the real non-ideal performance.


/***************************************************************************************************************/
Common:
[x] Add helper functions for parsing the input into Geometry, OperatingConditions, and ThermoProps structs

Impeller:
[x] Abstract Japikse function from ImpellerWindow class to Impeller class
[ ] Abstract the root finding to it's own math library
[ ] Add in Galvas model (from radcomp)
[ ] Add in CCMD model (in reference)

** bloackage is not being accounted for in GUI calculations even though it is in the formulas, need to double check
   what is going on there.