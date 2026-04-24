.. _download:

Installation
============

Download an executable file from the `executables <https://github.com/NewcastleRSE/RIA/tree/main/executables>`_ folder on Github for your system and off you go, or do the following.

|

1. Download the code from the home page.
2. Compile it by typing something like the following:

.. code-block:: none

    g++ -O3 src/*.cpp -o ria

3. It is also necessary to having working copies of the programs `PLINK <https://zzz.bwh.harvard.edu/plink/index.shtml>`_ and `KING version 2.2.9 <https://www.kingrelatedness.com/history.shtml>`_ **Note:** It is very important to download KING version 2.2.9 and not any later versions as we need to use the `--homog` option which is no longer available in later versions. 
4. Start analysing your data with RIA!
5. If you are using Windows then you will need to uncomment the line below in the `main.h` file:

.. code-block:: none

    #define USING_WINDOWS
