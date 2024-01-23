def turbine_calculation():
    while True:
        try:

            ###INPUT DATA###
            title = "The Content of Calculations"
            print(title, "\n", "-"*len(title), sep="")
            print("1-Turbine Power Calculation")
            print("2-TIV Selection Chart")
            print("3-Suction Head Calculation")
            print("4-Specific Speed Calculation")
            print("5-Runaway Speed Calculation")
            print("6-TIV Closing Time Calculation", "\n", "-"*len("6-TIV Closing Time Calculation"), sep="")
            print("\n")

            ##1-Turbine Power Calculation##
            print("1-Turbine Power Calculation", "\n", "-"*len("1-Turbine Power Calculation"), sep="")
            Q = float(input("Flow rate: [m3/s]..."))
            Hn = float(input("Net head: [m]..."))
            eta = float(input("Efficiency: [-]..."))
            rho = float(input("Density: [kg/m3]..."))
            g = float(input("Gravity: [m/s2]..."))

            Pturbine = rho*g*eta*Hn*Q / 1000
            #print(round(float(Pturbine)))

            print("Turbine Mechanical Power: {} [kW]".format(round(float(Pturbine))))
            print("\n")

            ##2-TIV Selection Chart##
            print("2-TIV Selection Chart", "\n", "-"*len("2-TIV Selection Chart"), sep="")

            def main():
                v = float(input("Reasonable Flow Velocity: [m/s]..."))

                while True:
                    Dpen = float(input("Penstock diameter: [m]..."))

                    Apen = 3.141592653589793 * (Dpen**2) / 4
                    print("Section Area in Penstock: {} [m^2]".format((float(Apen))))

                    Vpen = Q / Apen
                    print("Flow Velocity in Penstock: {} [m/s]".format((float(Vpen))))

                    if Vpen <= v:
                        print("Calculated flow velocity is acceptable")
                        break

                    else:
                        print("Warning: Calculated flow velocity in penstock is greater than reasonable flow velocity!")
                        to_continue = input("Press 'e' to restart the calculation, press 'h' to exit: ")

                        if to_continue.lower() != 'e':
                            break  # Exit the program!

            if __name__ == "__main__":
                main()

            print("\n")

            ##3-Suction Head Calculation##
            print("3-Suction Head Calculation", "\n", "-"*len("3-Suction Head Calculation"), sep="")

            Pmax = Pturbine #kW
            n = float(input("Revolution: [rpm]..."))
            D2 = float(input("Runner diameter: [m]..."))

            RunnerCL = float(input("Runner Centerline: [m]..."))
            TWL = float(input("Tail Water Level: [m]..."))

            hs = (RunnerCL - TWL)
            Hathd = 10.34 - 0.18 - 0.0012*RunnerCL

            sigma_plant = (Hathd - hs - (D2 / 2))/(Hn)
            sigma_turbine = 0.028

            safety = (sigma_plant / sigma_turbine)
            print("Safety factor: {} ".format((float(safety))))
            print("\n")

            ##4-Specific Speed Calculation##
            print("4-Specific Speed Calculation", "\n", "-"*len("4-Specific Speed Calculation"), sep="")

            nq_opt = n * (Q**(1/2))/ (Hn**(3/4))
            print("Specific speed in the optimum operating point: {} ".format(round(float(nq_opt))))

            ns_opt = 3.65 * nq_opt
            print("Specific speed related to the turbine performance: {} ".format(round(float(ns_opt))))
            print("\n")

            ##5-Runaway Speed Calculation##
            print("5-Runaway Speed Calculation", "\n", "-"*len("5-Runaway Speed Calculation"), sep="")
            psi = float(input("Minimum value at model runaway tests and varoius openings: psi..."))

            #Max steady state runaway speed
            nd= (((9.81*Hn) / (((3.141592653589793**2)/2)*(D2**2)*(psi)))**(1/2))*60
            print("Max steady state runaway speed: {} [rpm]".format(round(float(nd))))

            #Max steady state runaway speed incl. safety margin
            s = 1.05
            nds = nd * s
            print("Max steady state runaway speed incl. safety margin: {} [rpm]".format(round(float(nds))))
            print("\n")

            ##6-TIV Closing Time Calculation##
            print("6-TIV Closing Time Calculation", "\n", "-"*len("6-TIV Closing Time Calculation"), sep="")

            import math

            def calculate_Area(diameter):
                return (math.pi * diameter**2) / 4

            def calculate_L_A(diameter, length):
                return length / calculate_Area(diameter)

            def calculate_Water_Staring_Time(Q, Hn, L_A1):
                return L_A1 * Q / (9.81 * Hn)

            def calculate_Water_Velocity_at_Valve(Q, dvalve):
                return (4 * Q) / (dvalve**2 * math.pi)

            def calculate_Coefficient(c, Hn):
                return c / Hn**(1/2)

            def calculate_f1(coef):
                if 0.2 < coef < 0.4:
                    return ((coef - 0.2) / 0.2) * (9 - 18) + 18
                elif 0.4 < coef < 0.6:
                    return ((coef - 0.2) / 0.2) * (6 - 9) + 9
                elif 0.6 < coef < 0.8:
                    return ((coef - 0.6) / 0.2) * (4.5 - 6) + 6
                elif coef > 0.8:
                    return coef * 10

            def calculate_Max_Pressure(Hn, dhzul):
                return Hn * (1 + dhzul)

            def calculate_Required_Min_Linear_Closing_Time(f1, Tw, dhzul):
                return f1 * Tw / dhzul

            dvalve = float(input("Valve diameter [m]: "))

            def main():
                num_pipes = int(input("Enter the number of pipes: "))
                
                total_L_A = 0  # Initialize total_L_A to accumulate values for all pipes

                for i in range(1, num_pipes + 1):
                    print(f"\nPipe {i}:")
                    
                    d1 = float(input("Penstock_Diameter [m]: "))
                    L1 = float(input("Penstock_Length [m]: "))
                    
                    A1 = calculate_Area(d1)

                    L_A1 = calculate_L_A(d1, L1) 
                    total_L_A += L_A1  # Accumulate L_A for each pipe

                Tw = calculate_Water_Staring_Time(Q, Hn, total_L_A)

                c = calculate_Water_Velocity_at_Valve(Q, dvalve)
                coef = calculate_Coefficient(c, Hn)

                f1 = calculate_f1(coef)
                print("\n")
                dhzul = float(input("Admissible water hammer (e.g., 0.4 for 40%): "))
                Hmax = calculate_Max_Pressure(Hn, dhzul)
                Ts = calculate_Required_Min_Linear_Closing_Time(f1, Tw, dhzul)

                print("\nResults for TIV Closing Time Calculation:")
                print(f"Total L_A: {total_L_A} [1/m] (Sum of L_A for all pipes)")
                print(f"Water Staring Time, Tw: {Tw} [s]")
                print(f"Water Velocity at Valve, c: {c} [m/s]")
                print(f"Coefficient, coef: {coef} [m^0.5/s]")
                print(f"Max Pressure with adm. water hammer, Hmax: {Hmax} [m]")
                print(f"Required minimum linear closing time, Ts: {Ts} [s]")

            if __name__ == "__main__":
                main()


            break
            
        except ValueError:
            print("Incorrect entry! Please enter a numerical value!")
        except ZeroDivisionError:
            print("Incorrect entry! The second number cannot be zero!")
        except Exception as e:
            print("Something went wrong!", e)

turbine_calculation()

