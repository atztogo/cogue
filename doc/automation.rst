Architecture of automation
===========================

Cogue is composed of three objects, task, taskset, and controller
(``autocalc``). Task corresponds to a calculation. Taskset is a set of
tasks and tasksets. A series of tasksets can be contained in a
task. Complicated process is designed by nesting tasks and tasksets.
Controller handles tasks and tasksets with the help of a queueing system.

