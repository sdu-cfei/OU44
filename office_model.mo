within OU44;
package office_model
  model office_zone "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=1.0 "Internal wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";


    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
      "Indoor temperature setpoint [degC]"                         annotation (
        Placement(transformation(extent={{-240,-150},{-212,-122}}),
          iconTransformation(extent={{-234,-186},{-206,-158}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
              {234,-116}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-10,-48},{10,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=vent)
      annotation (Placement(transformation(extent={{-86,-76},{-66,-56}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMax=1,
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600)
      annotation (Placement(transformation(extent={{-150,-146},{-130,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput WS "VAV damper position [%]" annotation (
       Placement(transformation(extent={{-238,-80},{-210,-52}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-214},{236,-194}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Math.Gain scale1(k=100) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={160,-204})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-134,32})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{84,120},{
            90,120},{90,-60},{0,-60},{0,-48.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
            {46,-120},{46,120}},               color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix.port_b, air.port)
      annotation (Line(points={{10,-38},{46,-38},{46,120}}, color={191,0,0}));
    connect(Infltration.y, airMix.Vve) annotation (Line(points={{-65,-66},{-28,
            -66},{-28,-43},{-11,-43}}, color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
            120},{46,120}},
                       color={191,0,0}));
    connect(PID_temp.y, heating.u) annotation (Line(points={{-129,-136},{-129,
            -136},{-84,-136}}, color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-226,-136},{-190,-136},{-152,-136}},
                                                         color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-140,-180},{-140,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(PID_temp.y, scale1.u) annotation (Line(points={{-129,-136},{-106,-136},
            {-106,-204},{148,-204}}, color={0,0,127}));
    connect(scale1.y, vpos)
      annotation (Line(points={{171,-204},{226,-204}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,34},{-186,34},{-186,32},{-146,32}},
                                                     color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-123,32},{-68,32},{-68,34},{-12,34}},
                                                    color={0,0,127}));
    connect(WS, Infltration.u) annotation (Line(points={{-224,-66},{-88,-66}},
                             color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-34.6},{-10.8,-34.6}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
          Rectangle(
            extent={{-208,-104},{16,-160}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-208,-106},{-170,-116}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Heating")}),
      experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone;

  model office_zone_nowindow "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=1.0 "Internal wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
      "Indoor temperature setpoint [degC]"                         annotation (
        Placement(transformation(extent={{-240,-150},{-212,-122}}),
          iconTransformation(extent={{-234,-186},{-206,-158}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
              {234,-116}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-10,-48},{10,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,18},{-212,46}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=0)
      annotation (Placement(transformation(extent={{-86,-78},{-66,-58}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600,
      yMax=100)
      annotation (Placement(transformation(extent={{-174,-146},{-154,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput WS "VAV damper position [%]" annotation (
       Placement(transformation(extent={{-236,-82},{-208,-54}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-210},{236,-190}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-134,32})));
    Modelica.Blocks.Math.Gain scale(k=0.01)
      annotation (Placement(transformation(extent={{-126,-146},{-106,-126}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{84,120},{
            90,120},{90,-60},{0,-60},{0,-48.8}},                   color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
            {46,-120},{46,120}},               color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix.port_b, air.port)
      annotation (Line(points={{10,-38},{46,-38},{46,120}}, color={191,0,0}));
    connect(Infltration.y, airMix.Vve) annotation (Line(points={{-65,-68},{-28,
            -68},{-28,-43},{-11,-43}}, color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
            120},{46,120}},
                       color={191,0,0}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-226,-136},{-176,-136}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-164,-180},{-164,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,32},{-146,32}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-123,32},{-68,32},{-68,34},{-12,34}},
                                                    color={0,0,127}));
    connect(WS, Infltration.u) annotation (Line(points={{-222,-68},{-88,-68}},
                             color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-34.6},{-10.8,-34.6}}, color={0,0,127}));
    connect(PID_temp.y, scale.u)
      annotation (Line(points={{-153,-136},{-128,-136}}, color={0,0,127}));
    connect(scale.y, heating.u)
      annotation (Line(points={{-105,-136},{-84,-136}}, color={0,0,127}));
    connect(vpos, vpos)
      annotation (Line(points={{226,-200},{226,-200}}, color={0,0,127}));
    connect(PID_temp.y, vpos) annotation (Line(points={{-153,-136},{-146,-136},
            {-146,-200},{226,-200}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
          Rectangle(
            extent={{-208,-104},{16,-160}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-208,-106},{-170,-116}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Heating")}),
      experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_nowindow;

  model office_zone_MPC "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=1.0 "Internal wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{68,110},{88,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,
              -110},{234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=300 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-60,-48},{-40,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Blocks.Interfaces.RealInput Tsup
      "ventialtion supply air temperature"
      annotation (Placement(transformation(extent={{-246,-50},{-218,-22}}),
          iconTransformation(extent={{-240,-8},{-212,20}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2
      annotation (Placement(transformation(extent={{-122,-42},{-108,-28}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{68,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
            120},{46,120}},
                       color={191,0,0}));
    connect(Occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(Vair, airMix.Vve) annotation (Line(points={{-234,-70},{-80,-70},{
            -80,-43},{-61,-43}}, color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{-40,-38},{-22,-38}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{-2,-38},{46,-38},{46,120}}, color={191,0,0}));
    connect(Tsup, toKelvin2.Celsius) annotation (Line(points={{-232,-36},{-144,
            -36},{-144,-35},{-123.4,-35}}, color={0,0,127}));
    connect(toKelvin2.Kelvin, airMix.Tve) annotation (Line(points={{-107.3,-35},
            {-84.05,-35},{-84.05,-34.6},{-60.8,-34.6}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u) annotation (Line(points={{-12,
            -48},{-12,-58},{142,-58}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad) annotation (Line(points={{-12,-48},{
            -12,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{88,120},{180,120}}, color={0,0,127}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{88,120},{
            98,120},{98,-70},{-50,-70},{-50,-48.8}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_MPC;

  model office_zone_MPC_est_para "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-250,152},{-222,180}}),
          iconTransformation(extent={{-242,154},{-214,182}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-248,76},{-220,104}}),
          iconTransformation(extent={{-240,90},{-212,118}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=2.86 "External wall thermal resistance";
    parameter Real RInt=0.102 "Internal wall thermal resistance";
    parameter Real tmass=40 "Zone thermal mass factor [-]";
    parameter Real imass=199 "Zone internal thermal mass factor [-]";
    parameter Real shgc=2.82 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=1746 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 1746 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 140 "heat gains from occupancy [W]";

    parameter Real maxHeat=1310 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2
          /3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{68,110},{88,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{210,-68},{230,-48}}),   iconTransformation(extent={{214,-110},
              {234,-90}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=171.5 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            293.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix airMix
      annotation (Placement(transformation(extent={{-60,-48},{-40,-28}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1, initType=Modelica.Blocks.Types.Init.NoInit)
      annotation (Placement(transformation(extent={{144,-68},{164,-48}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-248,20},{-220,48}}),
          iconTransformation(extent={{-240,44},{-212,72}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=293.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(
          2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{210,-92},{230,-72}}),
          iconTransformation(extent={{216,-148},{236,-128}})));
    Modelica.Blocks.Interfaces.RealInput Vair "ventialtion air flow rate"
      annotation (Placement(transformation(extent={{-248,-84},{-220,-56}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-132,34})));
    Modelica.Blocks.Interfaces.RealInput Tsup
      "ventialtion supply air temperature"
      annotation (Placement(transformation(extent={{-246,-50},{-218,-22}}),
          iconTransformation(extent={{-240,-8},{-212,20}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor1
      annotation (Placement(transformation(extent={{-22,-48},{-2,-28}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin2
      annotation (Placement(transformation(extent={{-122,-42},{-108,-28}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-236,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-234,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{68,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{220,-58},{165,-58}},   color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,
            120},{46,120}},
                       color={191,0,0}));
    connect(Occ, gain.u)
      annotation (Line(points={{-234,34},{-144,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-121,34},{-12,34}}, color={0,0,127}));
    connect(Vair, airMix.Vve) annotation (Line(points={{-234,-70},{-80,-70},{-80,-43},
            {-61,-43}}, color={0,0,127}));
    connect(airMix.port_b, heatFlowSensor1.port_a)
      annotation (Line(points={{-40,-38},{-22,-38}}, color={191,0,0}));
    connect(heatFlowSensor1.port_b, air.port)
      annotation (Line(points={{-2,-38},{46,-38},{46,120}}, color={191,0,0}));
    connect(Tsup, toKelvin2.Celsius) annotation (Line(points={{-232,-36},{-144,-36},
            {-144,-35},{-123.4,-35}}, color={0,0,127}));
    connect(toKelvin2.Kelvin, airMix.Tve) annotation (Line(points={{-107.3,-35},{-84.05,
            -35},{-84.05,-34.6},{-60.8,-34.6}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, integrator.u)
      annotation (Line(points={{-12,-48},{-12,-58},{142,-58}}, color={0,0,127}));
    connect(heatFlowSensor1.Q_flow, qrad)
      annotation (Line(points={{-12,-48},{-12,-82},{220,-82}}, color={0,0,127}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{88,120},{180,120}}, color={0,0,127}));
    connect(temperatureSensor.T, airMix.Ti) annotation (Line(points={{88,120},{98,
            120},{98,-70},{-50,-70},{-50,-48.8}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-24},{16,-94}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-40},{-132,-56}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=18,
            textString="Ventialtion cooling
")}), experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_MPC_est_para;

  model office_zone_v2 "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=1.0 "Internal wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
      "Indoor temperature setpoint [degC]"                         annotation (
        Placement(transformation(extent={{-238,-150},{-210,-122}}),
          iconTransformation(extent={{-234,-186},{-206,-158}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
              {234,-116}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=45 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=vent)
      annotation (Placement(transformation(extent={{-124,-62},{-104,-42}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMax=1,
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600)
      annotation (Placement(transformation(extent={{-148,-146},{-128,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput Wineast "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-238,-66},{-210,-38}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-214},{236,-194}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Math.Gain scale1(k=100) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={160,-204})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-122,34})));
    Modelica.Blocks.Interfaces.RealInput Winwest "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-236,-94},{-208,-66}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Math.Gain Infltration1(k=vent)
      annotation (Placement(transformation(extent={{-124,-90},{-104,-70}})));
    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{-74,-72},{-54,-52}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{84,120},
            {90,120},{90,-60},{0,-60},{0,-50.8}}, color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
            {46,-120},{46,120}},               color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix_new.port_b, air.port)
      annotation (Line(points={{10,-40},{46,-40},{46,120}}, color={191,0,0}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{
            60,120},{46,120}},
                       color={191,0,0}));
    connect(PID_temp.y, heating.u) annotation (Line(points={{-127,-136},{-84,
            -136}},            color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-224,-136},{-150,-136}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-138,-180},{-138,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(PID_temp.y, scale1.u) annotation (Line(points={{-127,-136},{-106,
            -136},{-106,-204},{148,-204}},
                                     color={0,0,127}));
    connect(scale1.y, vpos)
      annotation (Line(points={{171,-204},{226,-204}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,34},{-134,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-111,34},{-12,34}}, color={0,0,127}));
    connect(Wineast, Infltration.u)
      annotation (Line(points={{-224,-52},{-126,-52}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix_new.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-36.6},{-10.8,-36.6}}, color={0,0,127}));
    connect(Winwest, Infltration1.u)
      annotation (Line(points={{-222,-80},{-126,-80}}, color={0,0,127}));
    connect(Infltration.y, add.u1) annotation (Line(points={{-103,-52},{-90,-52},{
            -90,-56},{-76,-56}}, color={0,0,127}));
    connect(Infltration1.y, add.u2) annotation (Line(points={{-103,-80},{-90,-80},
            {-90,-68},{-76,-68}}, color={0,0,127}));
    connect(add.y, airMix_new.Vve) annotation (Line(points={{-53,-62},{-32,-62},{-32,
            -45},{-11,-45}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,140},{-100,118}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)",
            textStyle={TextStyle.Bold}),
          Text(
            extent={{-202,56},{-88,42}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=14,
            textStyle={TextStyle.Bold},
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-206,218},{-132,196}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains",
            textStyle={TextStyle.Bold}),
          Rectangle(
            extent={{-208,-22},{16,-92}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-230,-12},{-100,-46}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            fontSize=14,
            textString="Window opening",
            textStyle={TextStyle.Bold}),
          Rectangle(
            extent={{-208,-104},{16,-160}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,-106},{-166,-116}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Heating",
            fontSize=14,
            textStyle={TextStyle.Bold})}),
      experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_v2;

  model office_zone_v2_R3C3 "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=1.0 "Internal wall thermal resistance";
    parameter Real Rwall=1.0 "zone wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real emass=10 "Zone wall thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 0 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re_wall(R=RExt/(3*Vi^
          (2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
      "Indoor temperature setpoint [degC]"                         annotation (
        Placement(transformation(extent={{-240,-150},{-212,-122}}),
          iconTransformation(extent={{-234,-186},{-206,-158}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
              {234,-116}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=45 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,122},{56,142}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=vent)
      annotation (Placement(transformation(extent={{-124,-62},{-104,-42}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMax=1,
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600)
      annotation (Placement(transformation(extent={{-150,-146},{-130,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput Wineast "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-238,-66},{-210,-38}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-214},{236,-194}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Math.Gain scale1(k=100) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={160,-204})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-122,34})));
    Modelica.Blocks.Interfaces.RealInput Winwest "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-236,-94},{-208,-66}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Math.Gain Infltration1(k=vent)
      annotation (Placement(transformation(extent={{-124,-90},{-104,-70}})));
    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{-74,-72},{-54,-52}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor Rintwall(R=Rwall/(3*
          Vi^(2/3)))
      annotation (Placement(transformation(extent={{-30,110},{-10,130}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air_wall(T(fixed=false,
          start=TAirInit + 273.15), C=emass*1.2*1005*Vi) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-76,84})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re_wall.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{56,120},{56,122},{46,122}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,122}},   color={191,0,0}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{84,120},
            {90,120},{90,-60},{0,-60},{0,-50.8}}, color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},{46,
            -120},{46,122}},                   color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix_new.port_b, air.port)
      annotation (Line(points={{10,-40},{46,-40},{46,122}}, color={191,0,0}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{60,122},
            {46,122}}, color={191,0,0}));
    connect(PID_temp.y, heating.u) annotation (Line(points={{-129,-136},{-129,
            -136},{-84,-136}}, color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-226,-136},{-190,-136},{-152,-136}},
                                                         color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-140,-180},{-140,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(PID_temp.y, scale1.u) annotation (Line(points={{-129,-136},{-106,-136},
            {-106,-204},{148,-204}}, color={0,0,127}));
    connect(scale1.y, vpos)
      annotation (Line(points={{171,-204},{226,-204}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,34},{-134,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-111,34},{-12,34}}, color={0,0,127}));
    connect(Wineast, Infltration.u)
      annotation (Line(points={{-224,-52},{-126,-52}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix_new.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-36.6},{-10.8,-36.6}}, color={0,0,127}));
    connect(Winwest, Infltration1.u)
      annotation (Line(points={{-222,-80},{-126,-80}}, color={0,0,127}));
    connect(Infltration.y, add.u1) annotation (Line(points={{-103,-52},{-90,-52},{
            -90,-56},{-76,-56}}, color={0,0,127}));
    connect(Infltration1.y, add.u2) annotation (Line(points={{-103,-80},{-90,-80},
            {-90,-68},{-76,-68}}, color={0,0,127}));
    connect(add.y, airMix_new.Vve) annotation (Line(points={{-53,-62},{-32,-62},{-32,
            -45},{-11,-45}}, color={0,0,127}));
    connect(re_wall.port_b, air_wall.port) annotation (Line(points={{-108,90},{
            -92,90},{-92,94},{-76,94}},
                                    color={191,0,0}));
    connect(air_wall.port, Rintwall.port_a)
      annotation (Line(points={{-76,94},{-76,120},{-30,120}}, color={191,0,0}));
    connect(Rintwall.port_b, air.port) annotation (Line(points={{-10,120},{18,120},
            {18,122},{46,122}}, color={191,0,0}));
    connect(prescribedHeatFlow.port, Rintwall.port_a) annotation (Line(points={{-138,
            166},{-58,166},{-58,120},{-30,120}}, color={191,0,0}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-22},{16,-92}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
          Rectangle(
            extent={{-208,-104},{16,-160}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-208,-106},{-170,-116}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Heating")}),
      experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_v2_R3C3;

  model office_zone_v2_with_para "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.351676819198899 "External wall thermal resistance";
    parameter Real RInt=0.03703283340918108 "Internal wall thermal resistance";
    parameter Real tmass=116.75499532735546 "Zone thermal mass factor [-]";
    parameter Real imass=395.11341897868544 "Zone internal thermal mass factor [-]";
    parameter Real shgc=2.296675041277071 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 170.6590043615286 "heat gains from occupancy [W]";

    parameter Real maxHeat=850.0068696767461 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-126,80},{-106,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput qh
      "Indoor temperature setpoint [degC]" annotation (Placement(transformation(
            extent={{-238,-146},{-210,-118}}), iconTransformation(extent={{-234,
              -186},{-206,-158}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
              {234,-116}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-80,-142},{-60,-122}})));
    parameter Real Vi=45 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            293.15))
      annotation (Placement(transformation(extent={{36,126},{56,146}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-122,34})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-60,-132},{-28,-132},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-126,90}}, color={191,0,0}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-106,90},{-80,90},{-80,126},{46,126}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,126},{46,126}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{56,120},{56,126},{46,126}},
                                                    color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,126}},   color={191,0,0}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
            {46,-120},{46,126}},               color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-134},{148,-134}},                color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{
            60,126},{46,126}},
                       color={191,0,0}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-172},{226,-172}},            color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,34},{-134,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-111,34},{-12,34}}, color={0,0,127}));
    connect(qh, prescribedHeatFlow1.Q_flow)
      annotation (Line(points={{-224,-132},{-80,-132}}, color={0,0,127}));
    connect(Tout, prescribedTemperature.T)
      annotation (Line(points={{-226,90},{-162,90}}, color={0,0,127}));
    connect(temperatureSensor.T, Tin)
      annotation (Line(points={{84,120},{224,120}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-22},{16,-92}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
          Rectangle(
            extent={{-208,-106},{16,-162}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-208,-106},{-170,-116}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Heating")}),
      experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_v2_with_para;

  model office_zone_v2_R1C1 "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=1.0 "Internal wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
      "Indoor temperature setpoint [degC]"                         annotation (
        Placement(transformation(extent={{-238,-150},{-210,-122}}),
          iconTransformation(extent={{-234,-186},{-206,-158}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
              {234,-116}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=45 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=vent)
      annotation (Placement(transformation(extent={{-124,-62},{-104,-42}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMax=1,
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600)
      annotation (Placement(transformation(extent={{-150,-146},{-130,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput Wineast "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-238,-66},{-210,-38}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-214},{236,-194}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Math.Gain scale1(k=100) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={160,-204})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-122,34})));
    Modelica.Blocks.Interfaces.RealInput Winwest "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-236,-94},{-208,-66}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Math.Gain Infltration1(k=vent)
      annotation (Placement(transformation(extent={{-124,-90},{-104,-70}})));
    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{-74,-72},{-54,-52}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{84,120},
            {90,120},{90,-60},{0,-60},{0,-50.8}}, color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
            {46,-120},{46,120}},               color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix_new.port_b, air.port)
      annotation (Line(points={{10,-40},{46,-40},{46,120}}, color={191,0,0}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(PID_temp.y, heating.u) annotation (Line(points={{-129,-136},{-129,
            -136},{-84,-136}}, color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-224,-136},{-152,-136}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-140,-180},{-140,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(PID_temp.y, scale1.u) annotation (Line(points={{-129,-136},{-106,-136},
            {-106,-204},{148,-204}}, color={0,0,127}));
    connect(scale1.y, vpos)
      annotation (Line(points={{171,-204},{226,-204}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,34},{-134,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-111,34},{-12,34}}, color={0,0,127}));
    connect(Wineast, Infltration.u)
      annotation (Line(points={{-224,-52},{-126,-52}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix_new.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-36.6},{-10.8,-36.6}}, color={0,0,127}));
    connect(Winwest, Infltration1.u)
      annotation (Line(points={{-222,-80},{-126,-80}}, color={0,0,127}));
    connect(Infltration.y, add.u1) annotation (Line(points={{-103,-52},{-90,-52},{
            -90,-56},{-76,-56}}, color={0,0,127}));
    connect(Infltration1.y, add.u2) annotation (Line(points={{-103,-80},{-90,-80},
            {-90,-68},{-76,-68}}, color={0,0,127}));
    connect(add.y, airMix_new.Vve) annotation (Line(points={{-53,-62},{-32,-62},{-32,
            -45},{-11,-45}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-22},{16,-92}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
          Rectangle(
            extent={{-208,-104},{16,-160}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-208,-106},{-170,-116}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Heating")}),
      experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_v2_R1C1;

  model office_zone_v2_R2C2_wall "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RInt=1.0 "Internal wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real shgc=0.5 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";

    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
      "Indoor temperature setpoint [degC]"                         annotation (
        Placement(transformation(extent={{-238,-150},{-210,-122}}),
          iconTransformation(extent={{-234,-186},{-206,-158}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
              {234,-116}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=45 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor wall(C=tmass*1.2*
          1005*Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{-74,88},{-54,108}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=vent)
      annotation (Placement(transformation(extent={{-124,-62},{-104,-42}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor zone(C=imass*1.2*
          1005*Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{32,120},{52,140}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-20,80},{0,100}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMax=1,
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600)
      annotation (Placement(transformation(extent={{-150,-146},{-130,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput Wineast "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-238,-66},{-210,-38}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-214},{236,-194}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Math.Gain scale1(k=100) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={160,-204})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-122,34})));
    Modelica.Blocks.Interfaces.RealInput Winwest "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-236,-94},{-208,-66}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Math.Gain Infltration1(k=vent)
      annotation (Placement(transformation(extent={{-124,-90},{-104,-70}})));
    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{-74,-72},{-54,-52}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{84,120},
            {90,120},{90,-60},{0,-60},{0,-50.8}}, color={0,0,127}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(re1.port_b, zone.port)
      annotation (Line(points={{0,90},{42,90},{42,120}}, color={191,0,0}));
    connect(PID_temp.y, heating.u) annotation (Line(points={{-129,-136},{-129,
            -136},{-84,-136}}, color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-224,-136},{-152,-136}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-140,-180},{-140,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(PID_temp.y, scale1.u) annotation (Line(points={{-129,-136},{-106,-136},
            {-106,-204},{148,-204}}, color={0,0,127}));
    connect(scale1.y, vpos)
      annotation (Line(points={{171,-204},{226,-204}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,34},{-134,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-111,34},{-12,34}}, color={0,0,127}));
    connect(Wineast, Infltration.u)
      annotation (Line(points={{-224,-52},{-126,-52}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix_new.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-36.6},{-10.8,-36.6}}, color={0,0,127}));
    connect(Winwest, Infltration1.u)
      annotation (Line(points={{-222,-80},{-126,-80}}, color={0,0,127}));
    connect(Infltration.y, add.u1) annotation (Line(points={{-103,-52},{-90,-52},{
            -90,-56},{-76,-56}}, color={0,0,127}));
    connect(Infltration1.y, add.u2) annotation (Line(points={{-103,-80},{-90,-80},
            {-90,-68},{-76,-68}}, color={0,0,127}));
    connect(add.y, airMix_new.Vve) annotation (Line(points={{-53,-62},{-32,-62},{-32,
            -45},{-11,-45}}, color={0,0,127}));
    connect(re.port_b, wall.port) annotation (Line(points={{-108,90},{-86,90},{
            -86,88},{-64,88}}, color={191,0,0}));
    connect(wall.port, re1.port_a) annotation (Line(points={{-64,88},{-42,88},{
            -42,90},{-20,90}}, color={191,0,0}));
    connect(zone.port, temperatureSensor.port)
      annotation (Line(points={{42,120},{64,120}}, color={191,0,0}));
    connect(prescribedHeatFlow.port, zone.port) annotation (Line(points={{-138,
            166},{-2,166},{-2,120},{42,120}}, color={191,0,0}));
    connect(occHeatGain.port, zone.port)
      annotation (Line(points={{8,34},{42,34},{42,120}}, color={191,0,0}));
    connect(airMix_new.port_b, zone.port)
      annotation (Line(points={{10,-40},{42,-40},{42,120}}, color={191,0,0}));
    connect(heatFlowSensor.port_b, zone.port) annotation (Line(points={{-4,-120},
            {42,-120},{42,120}}, color={191,0,0}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-22},{16,-92}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
          Rectangle(
            extent={{-208,-104},{16,-160}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-208,-106},{-170,-116}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Heating")}),
      experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_v2_R2C2_wall;

  model office_zone_v2_R5C4 "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));

    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.0 "External wall thermal resistance";
    parameter Real RExt_em=1.0 "External wall thermal resistance";
    parameter Real Rz=1.0 "External wall thermal resistance";
    parameter Real RInt=1.0 "Internal wall thermal resistance";
    parameter Real Rinf=1.0 "Internal wall thermal resistance";
    parameter Real tmass=5 "Zone thermal mass factor [-]";
    parameter Real imass=10 "Zone internal thermal mass factor [-]";
    parameter Real emass=10 "Zone internal thermal mass factor [-]";
    parameter Real emass_em=10 "Zone internal thermal mass factor [-]";
    parameter Real shgc1=0.5 "Solar heat gain coefficient [-]";
    parameter Real shgc2=0.5 "Solar heat gain coefficient [-]";

    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 150 "heat gains from occupancy [W]";
    parameter Real coeff1= 0.2 "heat gains from occupancy [W]";
    parameter Real coeff2= 0.4 "heat gains from occupancy [W]";
    parameter Real maxHeat=5000 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-130,80},{-110,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc2)
      annotation (Placement(transformation(extent={{-194,158},{-180,172}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain2
      annotation (Placement(transformation(extent={{-78,24},{-58,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
      "Indoor temperature setpoint [degC]"                         annotation (
        Placement(transformation(extent={{-238,-150},{-210,-122}}),
          iconTransformation(extent={{-234,-186},{-206,-158}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
              {234,-116}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=45 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{40,120},{60,140}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=vent)
      annotation (Placement(transformation(extent={{-124,-62},{-104,-42}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,160},{128,180}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor rin(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,150},{96,170}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMax=1,
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600)
      annotation (Placement(transformation(extent={{-150,-146},{-130,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput Wineast "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-238,-66},{-210,-38}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-214},{236,-194}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Math.Gain scale1(k=100) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={160,-204})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{6,-6},{-6,6}},
          rotation=180,
          origin={-154,34})));
    Modelica.Blocks.Interfaces.RealInput Winwest "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-236,-94},{-208,-66}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Math.Gain Infltration1(k=vent)
      annotation (Placement(transformation(extent={{-124,-90},{-104,-70}})));
    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{-74,-72},{-54,-52}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor ce(T(fixed=false,
          start=TAirInit + 273.15), C=emass*1.2*1005*Vi) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-84,100})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor rem(R=RExt_em/(3*Vi^(
          2/3)))
      annotation (Placement(transformation(extent={{-62,80},{-42,100}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor cem(T(fixed=false,
          start=TAirInit + 273.15), C=emass_em*1.2*1005*Vi) annotation (Placement(
          transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-28,100})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor rzone(R=Rz/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-12,80},{8,100}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor rinf(R=Rinf/(3*Vi^(2/
          3)))
      annotation (Placement(transformation(extent={{-66,110},{-46,130}})));
    Modelica.Blocks.Math.Gain solarCoeff1(k=shgc1)
      annotation (Placement(transformation(extent={{-194,184},{-180,198}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow2
      annotation (Placement(transformation(extent={{-158,182},{-138,202}})));
    Modelica.Blocks.Math.Gain factor(k=coeff1) annotation (Placement(
          transformation(
          extent={{5,-5},{-5,5}},
          rotation=180,
          origin={-107,47})));
    Modelica.Blocks.Math.Gain factor1(k=coeff2) annotation (Placement(
          transformation(
          extent={{5,-5},{-5,5}},
          rotation=180,
          origin={-107,31})));
    Modelica.Blocks.Math.Gain factor3(k=1 - coeff1 - coeff2) annotation (
        Placement(transformation(
          extent={{5,-5},{-5,5}},
          rotation=180,
          origin={-107,13})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain1
      annotation (Placement(transformation(extent={{-78,40},{-58,60}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain3
      annotation (Placement(transformation(extent={{-78,2},{-58,22}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-212,166},{
            -212,165},{-195.4,165}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-179.3,
            165},{-168,165},{-168,166},{-158,166}},
                                         color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-130,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{50,120}},  color={191,0,0}));
    connect(occHeatGain2.port, air.port)
      annotation (Line(points={{-58,34},{50,34},{50,120}}, color={191,0,0}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{84,120},
            {90,120},{90,-60},{0,-60},{0,-50.8}}, color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},{50,
            -120},{50,120}},                   color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix_new.port_b, air.port)
      annotation (Line(points={{10,-40},{50,-40},{50,120}}, color={191,0,0}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(rin.port_b, IntMass.port)
      annotation (Line(points={{96,160},{118,160}},           color={191,0,0}));
    connect(PID_temp.y, heating.u) annotation (Line(points={{-129,-136},{-129,
            -136},{-84,-136}}, color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-224,-136},{-152,-136}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-140,-180},{-140,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(PID_temp.y, scale1.u) annotation (Line(points={{-129,-136},{-106,-136},
            {-106,-204},{148,-204}}, color={0,0,127}));
    connect(scale1.y, vpos)
      annotation (Line(points={{171,-204},{226,-204}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,34},{-161.2,34}},
                                                     color={0,0,127}));
    connect(Wineast, Infltration.u)
      annotation (Line(points={{-224,-52},{-126,-52}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix_new.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-36.6},{-10.8,-36.6}}, color={0,0,127}));
    connect(Winwest, Infltration1.u)
      annotation (Line(points={{-222,-80},{-126,-80}}, color={0,0,127}));
    connect(Infltration.y, add.u1) annotation (Line(points={{-103,-52},{-90,-52},{
            -90,-56},{-76,-56}}, color={0,0,127}));
    connect(Infltration1.y, add.u2) annotation (Line(points={{-103,-80},{-90,-80},
            {-90,-68},{-76,-68}}, color={0,0,127}));
    connect(add.y, airMix_new.Vve) annotation (Line(points={{-53,-62},{-32,-62},{-32,
            -45},{-11,-45}}, color={0,0,127}));
    connect(re.port_b, ce.port)
      annotation (Line(points={{-110,90},{-84,90}}, color={191,0,0}));
    connect(ce.port, rem.port_a)
      annotation (Line(points={{-84,90},{-62,90}}, color={191,0,0}));
    connect(rem.port_b, cem.port)
      annotation (Line(points={{-42,90},{-28,90}}, color={191,0,0}));
    connect(cem.port, rzone.port_a)
      annotation (Line(points={{-28,90},{-12,90}}, color={191,0,0}));
    connect(rzone.port_b, air.port) annotation (Line(points={{8,90},{28,90},{28,120},
            {50,120}}, color={191,0,0}));
    connect(rin.port_a, air.port) annotation (Line(points={{76,160},{64,160},{64,120},
            {50,120}}, color={191,0,0}));
    connect(prescribedTemperature.port, rinf.port_a) annotation (Line(points={{-140,
            90},{-134,90},{-134,120},{-66,120}}, color={191,0,0}));
    connect(rinf.port_b, air.port)
      annotation (Line(points={{-46,120},{50,120}}, color={191,0,0}));
    connect(prescribedHeatFlow.port, ce.port) annotation (Line(points={{-138,166},
            {-96,166},{-96,78},{-84,78},{-84,90}}, color={191,0,0}));
    connect(Solrad, solarCoeff1.u) annotation (Line(points={{-228,166},{-214,166},
            {-214,191},{-195.4,191}}, color={0,0,127}));
    connect(solarCoeff1.y, prescribedHeatFlow2.Q_flow) annotation (Line(points={{-179.3,
            191},{-168.65,191},{-168.65,192},{-158,192}}, color={0,0,127}));
    connect(prescribedHeatFlow2.port, IntMass.port) annotation (Line(points={{-138,
            192},{104,192},{104,160},{118,160}}, color={191,0,0}));
    connect(factor.y, occHeatGain1.Q_flow) annotation (Line(points={{-101.5,47},{-90.75,
            47},{-90.75,50},{-78,50}}, color={0,0,127}));
    connect(factor1.y, occHeatGain2.Q_flow) annotation (Line(points={{-101.5,31},{
            -89.75,31},{-89.75,34},{-78,34}}, color={0,0,127}));
    connect(factor3.y, occHeatGain3.Q_flow) annotation (Line(points={{-101.5,13},{
            -90.75,13},{-90.75,12},{-78,12}}, color={0,0,127}));
    connect(gain.y, factor.u) annotation (Line(points={{-147.4,34},{-130,34},{-130,
            47},{-113,47}}, color={0,0,127}));
    connect(gain.y, factor1.u) annotation (Line(points={{-147.4,34},{-130,34},{-130,
            31},{-113,31}}, color={0,0,127}));
    connect(gain.y, factor3.u) annotation (Line(points={{-147.4,34},{-130,34},{-130,
            13},{-113,13}}, color={0,0,127}));
    connect(occHeatGain1.port, cem.port)
      annotation (Line(points={{-58,50},{-28,50},{-28,90}}, color={191,0,0}));
    connect(occHeatGain3.port, IntMass.port)
      annotation (Line(points={{-58,12},{118,12},{118,160}}, color={191,0,0}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-22},{16,-92}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
          Rectangle(
            extent={{-208,-104},{16,-160}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-208,-106},{-170,-116}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Heating")}),
      experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_v2_R5C4;

  model office_zone_v2_rbc "Simple zone"

    package Air = Buildings.Media.Air(extraPropertiesNames={"CO2"})
      "Medium model";

    Modelica.Blocks.Interfaces.RealInput Solrad "Solar radiation [W]"
      annotation (Placement(transformation(extent={{-242,152},{-214,180}}),
          iconTransformation(extent={{-234,166},{-206,194}})));
    Modelica.Blocks.Interfaces.RealInput Tout "Outdoor temperature [degC]"
      annotation (Placement(transformation(extent={{-240,76},{-212,104}}),
          iconTransformation(extent={{-234,84},{-206,112}})));
    Modelica.Blocks.Interfaces.RealOutput Tin " degCelsius" annotation (
        Placement(transformation(extent={{214,110},{234,130}}),
          iconTransformation(extent={{214,110},{234,130}})));
    parameter Real TAirInit=20 "Initial temperature of indoor air [degC]";
    parameter Real CO2Init=400 "Initial CO2 concentration [ppmv]";
    parameter Real RExt=1.9027346543788286 "External wall thermal resistance";
    parameter Real RInt=0.05920884641183984 "Internal wall thermal resistance";
    parameter Real tmass=199.08618970594372 "Zone thermal mass factor [-]";
    parameter Real imass=952.9802343295498 "Zone internal thermal mass factor [-]";
    parameter Real shgc=2.2218134182575815 "Solar heat gain coefficient [-]";
    parameter Real CO2n=400 "CO2 Neutral Level";
    parameter Real CO2pp(unit="m3/h")=0.02 "CO2 generation per person [m3/h]";
    parameter Real maxVent=2000 "Maximum ventilation flowrate [m3/h]";
    parameter Real vent= 2000 "air infltration from window opening[m3/h]";
    parameter Real occ_gain= 196.36645064908382 "heat gains from occupancy [W]";

    parameter Real maxHeat=850.0068696767461 "Heating power of radiators [W]";

    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re(R=RExt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{-128,80},{-108,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
      annotation (Placement(transformation(extent={{64,110},{84,130}})));
    Modelica.Blocks.Math.Gain solarCoeff(k=shgc)
      annotation (Placement(transformation(extent={{-196,156},{-176,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
      annotation (Placement(transformation(extent={{-158,156},{-138,176}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow occHeatGain
      annotation (Placement(transformation(extent={{-12,24},{8,44}})));
    Modelica.Blocks.Interfaces.RealInput Tset
      "Indoor temperature setpoint [degC]"                         annotation (
        Placement(transformation(extent={{-238,-150},{-210,-122}}),
          iconTransformation(extent={{-234,-186},{-206,-158}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature
      prescribedTemperature
      annotation (Placement(transformation(extent={{-160,80},{-140,100}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heatFlowSensor
      annotation (Placement(transformation(extent={{-24,-130},{-4,-110}})));
    Modelica.Blocks.Interfaces.RealOutput Qrad
      "Heat supplied by radiators [Wh]"
                                       annotation (Placement(transformation(
            extent={{216,-144},{236,-124}}), iconTransformation(extent={{214,-136},
              {234,-116}})));
    Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
      annotation (Placement(transformation(extent={{-54,-146},{-34,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin
      annotation (Placement(transformation(extent={{182,110},{202,130}})));
    Modelica.Thermal.HeatTransfer.Celsius.ToKelvin toKelvin1
      annotation (Placement(transformation(extent={{-200,80},{-180,100}})));
    parameter Real Vi=45 "Air volume [m3]";
    parameter Real occheff=1. "Occupant heat generation effectiveness";
    parameter Real Vinf=50 "Infiltration rate";
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor air(
                                    C=tmass*1.2*1005*Vi, T(fixed=false, start=
            TAirInit + 273.15))
      annotation (Placement(transformation(extent={{36,120},{56,140}})));
    Components.AirMix_new
                      airMix_new
      annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
    Modelica.Blocks.Continuous.Integrator integrator(k=1)
      annotation (Placement(transformation(extent={{150,-144},{170,-124}})));
    Modelica.Blocks.Interfaces.RealInput Occ "Occupancy existance"
      annotation (Placement(transformation(extent={{-240,20},{-212,48}}),
          iconTransformation(extent={{-234,2},{-206,30}})));
    Modelica.Blocks.Math.Gain Infltration(k=vent)
      annotation (Placement(transformation(extent={{-124,-62},{-104,-42}})));
    Modelica.Blocks.Math.Gain heating(k=maxHeat)
      annotation (Placement(transformation(extent={{-82,-146},{-62,-126}})));
    Modelica.Thermal.HeatTransfer.Components.HeatCapacitor IntMass(C=imass*1.2*1005
          *Vi, T(fixed=false, start=TAirInit + 273.15))
      annotation (Placement(transformation(extent={{108,156},{128,176}})));
    Modelica.Thermal.HeatTransfer.Components.ThermalResistor re1(R=RInt/(3*Vi^(2/3)))
      annotation (Placement(transformation(extent={{76,148},{96,168}})));
    Modelica.Blocks.Continuous.LimPID PID_temp(
      yMax=1,
      yMin=0,
      initType=Modelica.Blocks.Types.InitPID.InitialState,
      k=0.5,
      controllerType=Modelica.Blocks.Types.SimpleController.PI,
      Ti=600)
      annotation (Placement(transformation(extent={{-148,-146},{-128,-126}})));
    Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin1
      annotation (Placement(transformation(extent={{-50,-190},{-70,-170}})));
    Modelica.Blocks.Interfaces.RealOutput qrad "Heat supply by radiators [W]"
      annotation (Placement(transformation(extent={{216,-182},{236,-162}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Interfaces.RealInput Wineast "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-238,-66},{-210,-38}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Interfaces.RealOutput vpos "Radiator valve position[%]"
      annotation (Placement(transformation(extent={{216,-214},{236,-194}}),
          iconTransformation(extent={{214,-136},{234,-116}})));
    Modelica.Blocks.Math.Gain scale1(k=100) annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={160,-204})));
    parameter Real Tve=21 "Ventilation air temperature";
    Modelica.Blocks.Math.Gain gain(k=occ_gain) annotation (Placement(
          transformation(
          extent={{10,-10},{-10,10}},
          rotation=180,
          origin={-122,34})));
    Modelica.Blocks.Interfaces.RealInput Winwest "VAV damper position [%]"
      annotation (Placement(transformation(extent={{-236,-94},{-208,-66}}),
          iconTransformation(extent={{-240,-60},{-212,-32}})));
    Modelica.Blocks.Math.Gain Infltration1(k=vent)
      annotation (Placement(transformation(extent={{-124,-90},{-104,-70}})));
    Modelica.Blocks.Math.Add add
      annotation (Placement(transformation(extent={{-74,-72},{-54,-52}})));
  equation
    connect(Solrad, solarCoeff.u) annotation (Line(points={{-228,166},{-198,166}},
                        color={0,0,127}));
    connect(solarCoeff.y, prescribedHeatFlow.Q_flow) annotation (Line(points={{-175,
            166},{-158,166}},            color={0,0,127}));
    connect(prescribedHeatFlow1.port, heatFlowSensor.port_a) annotation (Line(
          points={{-34,-136},{-28,-136},{-28,-120},{-24,-120}},  color={191,0,0}));
    connect(prescribedTemperature.port, re.port_a)
      annotation (Line(points={{-140,90},{-128,90}}, color={191,0,0}));
    connect(temperatureSensor.T, fromKelvin.Kelvin)
      annotation (Line(points={{84,120},{180,120}},            color={0,0,127}));
    connect(Tin, fromKelvin.Celsius) annotation (Line(points={{224,120},{210,
            120},{203,120}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, prescribedTemperature.T)
      annotation (Line(points={{-179,90},{-162,90}},           color={0,0,127}));
    connect(toKelvin1.Celsius,Tout)
      annotation (Line(points={{-202,90},{-226,90}},           color={0,0,127}));
    connect(re.port_b, air.port)
      annotation (Line(points={{-108,90},{-80,90},{-80,120},{46,120}},
                                                   color={191,0,0}));
    connect(prescribedHeatFlow.port, air.port) annotation (Line(points={{-138,
            166},{-80,166},{-80,120},{46,120}},
                                           color={191,0,0}));
    connect(temperatureSensor.port, air.port)
      annotation (Line(points={{64,120},{46,120}},  color={191,0,0}));
    connect(occHeatGain.port, air.port)
      annotation (Line(points={{8,34},{46,34},{46,120}},   color={191,0,0}));
    connect(temperatureSensor.T, airMix_new.Ti) annotation (Line(points={{84,120},
            {90,120},{90,-60},{0,-60},{0,-50.8}}, color={0,0,127}));
    connect(heatFlowSensor.port_b, air.port) annotation (Line(points={{-4,-120},
            {46,-120},{46,120}},               color={191,0,0}));
    connect(heatFlowSensor.Q_flow, integrator.u) annotation (Line(points={{-14,
            -130},{-14,-130},{-14,-134},{148,-134}},     color={0,0,127}));
    connect(Qrad, integrator.y)
      annotation (Line(points={{226,-134},{176,-134},{171,-134}},
                                                       color={0,0,127}));
    connect(airMix_new.port_b, air.port)
      annotation (Line(points={{10,-40},{46,-40},{46,120}}, color={191,0,0}));
    connect(prescribedHeatFlow1.Q_flow, heating.y) annotation (Line(points={{-54,
            -136},{-54,-136},{-61,-136}},        color={0,0,127}));
    connect(re1.port_b, IntMass.port)
      annotation (Line(points={{96,158},{118,158},{118,156}}, color={191,0,0}));
    connect(re1.port_a, air.port) annotation (Line(points={{76,158},{60,158},{
            60,120},{46,120}},
                       color={191,0,0}));
    connect(PID_temp.y, heating.u) annotation (Line(points={{-127,-136},{-84,
            -136}},            color={0,0,127}));
    connect(Tset, PID_temp.u_s)
      annotation (Line(points={{-224,-136},{-150,-136}}, color={0,0,127}));
    connect(fromKelvin1.Celsius, PID_temp.u_m) annotation (Line(points={{-71,
            -180},{-138,-180},{-138,-148}}, color={0,0,127}));
    connect(fromKelvin1.Kelvin, temperatureSensor.T) annotation (Line(points={{-48,
            -180},{90,-180},{90,120},{84,120}},     color={0,0,127}));
    connect(heatFlowSensor.Q_flow, qrad) annotation (Line(points={{-14,-130},{
            -14,-130},{-14,-172},{226,-172}}, color={0,0,127}));
    connect(PID_temp.y, scale1.u) annotation (Line(points={{-127,-136},{-106,
            -136},{-106,-204},{148,-204}},
                                     color={0,0,127}));
    connect(scale1.y, vpos)
      annotation (Line(points={{171,-204},{226,-204}}, color={0,0,127}));
    connect(Occ, gain.u)
      annotation (Line(points={{-226,34},{-134,34}}, color={0,0,127}));
    connect(gain.y, occHeatGain.Q_flow)
      annotation (Line(points={{-111,34},{-12,34}}, color={0,0,127}));
    connect(Wineast, Infltration.u)
      annotation (Line(points={{-224,-52},{-126,-52}}, color={0,0,127}));
    connect(toKelvin1.Kelvin, airMix_new.Tve) annotation (Line(points={{-179,90},{
            -176,90},{-176,-36.6},{-10.8,-36.6}}, color={0,0,127}));
    connect(Winwest, Infltration1.u)
      annotation (Line(points={{-222,-80},{-126,-80}}, color={0,0,127}));
    connect(Infltration.y, add.u1) annotation (Line(points={{-103,-52},{-90,-52},{
            -90,-56},{-76,-56}}, color={0,0,127}));
    connect(Infltration1.y, add.u2) annotation (Line(points={{-103,-80},{-90,-80},
            {-90,-68},{-76,-68}}, color={0,0,127}));
    connect(add.y, airMix_new.Vve) annotation (Line(points={{-53,-62},{-32,-62},{-32,
            -45},{-11,-45}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-220,
              -220},{220,220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-208,138},{-100,66}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-208,62},{16,-18}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-202,136},{-104,122}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="External walls (cond. and inf.)"),
          Text(
            extent={{-204,54},{-118,40}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            fontSize=18,
            textString="Occ. heat gains
"),       Rectangle(
            extent={{-208,214},{-124,146}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-204,212},{-150,202}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            horizontalAlignment=TextAlignment.Left,
            textString="Solar heat gains"),
          Rectangle(
            extent={{-208,-22},{16,-92}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-210,-38},{-132,-54}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Window opening",
            fontSize=18),
          Rectangle(
            extent={{-208,-104},{16,-160}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-208,-106},{-170,-116}},
            lineColor={28,108,200},
            fillColor={225,225,225},
            fillPattern=FillPattern.Solid,
            textString="Heating")}),
      experiment(Tolerance=1e-006),
      __Dymola_experimentSetupOutput,
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-220,-220},{220,
              220}},
          initialScale=0.1), graphics={
          Rectangle(
            extent={{-220,-220},{220,220}},
            lineColor={95,95,95},
            fillColor={95,95,95},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-202,202},{200,-202}},
            pattern=LinePattern.None,
            lineColor={117,148,176},
            fillColor={170,213,255},
            fillPattern=FillPattern.Sphere),
          Rectangle(
            extent={{-96,102},{96,-100}},
            lineColor={0,0,0},
            fillColor={135,135,135},
            fillPattern=FillPattern.Solid)}));
  end office_zone_v2_rbc;
end office_model;
