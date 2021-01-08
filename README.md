# 8021aqView
8021aq-View is a text-based python tool, to simulate and understand Shortest Path Bridging networks and analyze path usage,
multicast states, computation time and much more.

![Example of the script](/images/screenshot1.png)

The tool is based on the files “8021aqView-v1.py” and “spbnetworks.py”, that needs to be placed in the same folder. In the second mentioned file, you can find three example networks – a super spine, the diamond used in RFC 6329 and a simple square topology. You can find the drawings in the folder "images".

You can create your own networks, based on those examples. Just create your nodes like in the given examples and adapt name, Bridge ID, Links and ISID´s.
The format for the key “Links” is:

<b>[[“Name of adjacent node”, “Metric”]]</b>

The format for the key “ISID” is:

<b>[[“ISID#”, “BVLAN”, “Multicast transmit bit”, “Multicast receive bit”]]</b>

Hence, if you would like your ISID to participate in multicast tandem S,G mode, you need to put both bits to 1. For multicast headend replication, set the bits to 0.

![Example of the script](/images/screenshot2.png)

You can create several different SPB networks in ”spbnetworks.py”, but only one is imported by “8021aqView-v1.py” and computed for further analysis. The network/dictionary you import, can be changed under in the bottom line of “spbnetworks.py”. Just uncomment the network you want to load:

![Example of the script](/images/screenshot3.png)

Start “8021awView-v1.py” and try the different options in the main menu, and change your SPB networks i.e. link metric, BVLANs, BridgeID etc. to see how it behaves. After every change in “spbnetworks.py”, the script needs to be restarted to take the changes into account. Don´t forget to save your changes in “spbnetworks.py”.
You can also set the “SystemState” in the dictionary to “off”. This will compute the SPB network as if the node would be physically down.
If your network contains a faulty configuration like duplicated BridgeID, ISID/BVLAN mismatch, unidirectional links etc., an exception will stop the start of the script and tell you what needs to be fixed.
You can also add additional “random” ISID`s/BVLAN`s by modifying the dictionary in the “8021aqView-v1.py”. This way, you can i.e. check how your computation time and multicast tables are impacted. Just uncomment the lines as needed, another 1.000 ISID`s (500 in headend mode and 500 in tandem mode) spread over 256 BVLAN`s are added in the example below.

![Example of the script](/images/screenshot4.png)

Be patient, if you add a lot of tandem ISID´s. Depending on your system, it can take a remarkable time to compute them and populate the multicast table. 
Even if the standard 802.1aq offers only 16 ECT values, to tool makes an experimental use of all 255 possible ECT-masks – to analyze maximal theoretical path diversity.

Please get in touch with me, if you have any questions or enhancement requests.




