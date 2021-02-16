# -*- mode: Make -*-

test:
	make test-refdata
	make test-functionality


test-stable:
	make test-refdata
	make test-functionality-stable

speed:
	cd speed && make speed_check


######################
test-refdata:
	cd refdata && make test

test-functionality:
	cd e2e && make test

test-functionality-stable:
	cd e2e && make test-stable

