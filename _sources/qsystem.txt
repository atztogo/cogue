.. _qsystem:

Installation and configure of queueing system
==============================================

Queueing system is one of the most important component of the Cogue
package. Probably researchers around first-principles calculation are
familiar with how to use queueing systems, however may not have
experience how to build the queueing system. In fact, it is not very
difficult to install and configure a queueing system. In this section,
installation and configuration of the Grid Engine software are
explained to take Ubuntu linux as an example of the operating system.


Installation of Grid Engine
-----------------------------

Grid Engine is in the list of Ubuntu packages. So the installation is
very simple. On the server node::

   % sudo apt-get install gridengine-client gridengine-common \
     gridengine-master gridengine-qmon

On the client nodes::

   % sudo apt-get install gridengine-exec gridengine-client

During the installation, you are asked several questions. For postfix,
you can choose no configuration. For Grid Engine, 'Configure Grid
Engine automatically' is Yes. Cell name is arbitrary and 'default' is
OK. Master hostname has to be resolved by ``/etc/hosts`` or DNS. A few
initial settings are stored under ``/var/lib/gridengine`` in the case
of Ubuntu linux.

If there is only one computer, it can behave as both server and
client simultaneously.

Configuration of Grid Engine
-----------------------------

Grid Engine has a graphical user interface for the configuration
called 'qmon'. On the master node, it is opened by::

   % sudo qmon

If the master node is located on the remote place, ``ssh -X`` is used
to login the master node. Then the graphics will be forwarded to your
computer. After opening qmon, if qmon reports some
authorization problem, then you should check the hostname written in
``/var/lib/gridengine/CELL_NAME/common/act_qmaster``. If you modify
the hostname in it, you need to restart gridendinge-master by::

   % sudo service gridengine-master restart


In most cases, to add to something already exists can be done from
"Modify" button, not from "Add" in qmon.  The first steps of the
configurations are as follows:

1. In "Host Configuration" -> "Execution Host" tab, execution hosts
   are added. The master host is probably already set as an
   administration host and a submit host.

2. In "User Configuration" -> "Userset" tab, users are added to the
   userset "arusers".

3. In "Queue Control" -> "Cluster Queues" tab, a new cluster queue is
   added.

   * In "General Configuration" tab, the number of cores in an
     execution node is specified in "Slots" field.
   
   * In "Use Access" tab, "arusers" is added to "Allow Access to"
     field.

   * If you have already created a PE (parallel environment) object,
     in "Parallel Environment" tab, the PE may be added to
     "Referenced PEs" field. The configuration of PE is tricky and
     will be explained in the next section.

The detailed information is found `here
<http://docs.oracle.com/cd/E24901_01/index.htm>`_.
`This web site <http://pka.engr.ccny.cuny.edu/~jmao/node/49>`_ is 
helpful when starting configuration for the first time.


Configuration of parallel environment
---------------------------------------

If you want to use multi-cores on CPU or multiple computers in MPI,
you have to configure the parallel environment of grid engine.
