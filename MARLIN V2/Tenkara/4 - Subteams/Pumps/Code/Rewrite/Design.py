from Pump import Pump

def Design(pump: Pump):
	import Pump1D
	print(Pump1D.GuessInducerlessNPSH_r(pump))
	print(pump.NPSH_a)

	inducer_head = Pump1D.GetInducerHeadRequirement(pump)
	print(inducer_head)
	pump.Impeller.DP.H_rise = pump.DesignHeadRise - inducer_head
	# TODO create inducer

	pump.Impeller.PrintDimensionlessParameters()


if __name__ == "__main__":
	from Main import LOXPump
	Design(LOXPump())