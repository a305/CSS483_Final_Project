package com.uwb.bt2j.aligner;

import com.uwb.bt2j.util.EList;

import javafx.util.Pair;

public class Presets {
	public Presets() {
		
	}
	
	public void apply(String preset, String policy, EList<Pair<Integer, String>> opts) {
		// Presets:                 Same as:
		//  For --end-to-end:
		//   --very-fast            -M 5 -R 1 -N 0 -L 22 -i S,1,2.50
		//   --fast                 -M 10 -R 2 -N 0 -L 22 -i S,1,2.50
		//   --sensitive            -M 15 -R 2 -N 0 -L 22 -i S,1,1.15
		//   --very-sensitive       -M 25 -R 3 -N 0 -L 19 -i S,1,0.50
		if(preset == "very-fast") {
			policy += ";SEED=0";
			policy += ";SEEDLEN=22";
			policy += ";DPS=5";
			policy += ";ROUNDS=1";
			policy += ";IVAL=S,0,2.50";
		} else if(preset == "fast") {
			policy += ";SEED=0";
			policy += ";SEEDLEN=22";
			policy += ";DPS=10";
			policy += ";ROUNDS=2";
			policy += ";IVAL=S,0,2.50";
		} else if(preset == "sensitive") {
			policy += ";SEED=0";
			policy += ";SEEDLEN=22";
			policy += ";DPS=15";
			policy += ";ROUNDS=2";
			policy += ";IVAL=S,1,1.15";
		} else if(preset == "very-sensitive") {
			policy += ";SEED=0";
			policy += ";SEEDLEN=20";
			policy += ";DPS=20";
			policy += ";ROUNDS=3";
			policy += ";IVAL=S,1,0.50";
		}
		//  For --local:
		//   --very-fast-local      -M 1 -N 0 -L 25 -i S,1,2.00
		//   --fast-local           -M 2 -N 0 -L 22 -i S,1,1.75
		//   --sensitive-local      -M 2 -N 0 -L 20 -i S,1,0.75 (default)
		//   --very-sensitive-local -M 3 -N 0 -L 20 -i S,1,0.50
		else if(preset == "very-fast-local") {
			policy += ";SEED=0";
			policy += ";SEEDLEN=25";
			policy += ";DPS=5";
			policy += ";ROUNDS=1";
			policy += ";IVAL=S,1,2.00";
		} else if(preset == "fast-local") {
			policy += ";SEED=0,22";
			policy += ";DPS=10";
			policy += ";ROUNDS=2";
			policy += ";IVAL=S,1,1.75";
		} else if(preset == "sensitive-local") {
			policy += ";SEED=0";
			policy += ";SEEDLEN=20";
			policy += ";DPS=15";
			policy += ";ROUNDS=2";
			policy += ";IVAL=S,1,0.75";
		} else if(preset == "very-sensitive-local") {
			policy += ";SEED=0";
			policy += ";SEEDLEN=20";
			policy += ";DPS=20";
			policy += ";ROUNDS=3";
			policy += ";IVAL=S,1,0.50";
		}
		else {
			System.err.println("Unknown preset: " + preset + "\n");
		}
	}
}
