pyforwind
=========

An open-source package to generate synthetic IEC-conform wind fields with extended turbulence characteristics. 

Installation
------------

The ``pyforwind`` package is available on pypi and can be installed using pip

.. code-block:: shell

    pip install pyforwind

How to use this package
-----------------------

Generate superstatistical Kaimal wind field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To use ``pyforwind`` to generate a Kaimal wind field ``u_kaimal`` and its superstatistical (i.e., non-Gaussian)
extension ``u_super_kaimal``, import the ``SFW`` with the model parameters: integral length scale ``L``, intermittency coefficient ``mu``,
horizontal wind speed at rotor hub ``V_hub``, hub height ``h_hub``, time length and diameter ``(T, diam)``, resolution ``(N_T, N_rotor)``,
and the wind field type ``kind``.

.. code-block:: python

    from pyforwind import SFW

    swf_kaimal = SWF(L, mu, V_hub, h_hub, (T, diam), (N_T, N_rotor), kind='gauss')
    swf_super_kaimal = SWF(L, mu, V_hub, h_hub, (T, diam), (N_T, N_rotor), kind='spatiotemporal')
    u_kaimal = swf_kaimal.field(seed)
    u_super_kaimal = swf_super_kaimal.field(seed)
