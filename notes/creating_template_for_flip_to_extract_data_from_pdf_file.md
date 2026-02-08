# Creating a Template for the Flip program to Extract Data from a PDF File

This file will walk you through creating a template from scratch for Flip.

## Switch off Multiline Mode on regex101.com

There is an important tip you need to know when making a new template for a PDF file in the Flip program. When you are using regex101.com to test and build the regex with your sample data, you need to switch of "multiline mode". By default, regex101 runs in "global" and "multiline mode" (flags gm). Turn off m flag and it will match how the C++ regex engine in Flip works. If you don't, you will get problems when you use your regexes in Flip and it will confuse you.

## Get your source data from Flip and put it into regex101.com

Open Flip, load a PDF file into it, tick "Show data before regex processing" and "Strip excessive whitespace", select any template at first (it doesn't matter you will make a new one for whatever data you're trying to get form the PDF) then click LAUNCH button. In the new window that appears click in the left panel where the PDF text data is, pres Ctrl+A to select all the text, Ctrl+C to copy, then paste that into the regex101.com input sample area. Now you are working on the exact text that the Flip program is working on and you can develop your regexes from here.

## Set regex101.com to "Substitution" Mode

Flip uses regex_replace or "substitution" as it's called on regex101.com. To put regex101.com into substitution mode, towards the left side of the page, you will see an area that says "Save & Share", "Flavor", "Function". In the Function section, click on "Substitution" and the middle panel will split to have a second pane at hte bottom and an extra input field. This is where you can use regex substitution tokens to alter your PDF data with the regex matching system. Sub tokens are like: $1, $2, $3 etc depending on what matches your regex does. This document does not go into details about how to do this. Look up a youtube about it.

## Create a new template in Flip

Open the Template Manager in Flip by focusing the main Flip window and pressing Ctrl+T. In the Template Manager window you can create a new template by just entering a new template filename and clicking the Add button. You can browse to a specific location or just type the name in and it will save in your user's home directory in the .flip folder (where all Flip program data for your user will be saved for easy backup)

You can also choose an existing template and edit it if that's what you'd prefer.

Templates save automatically a few seconds after you stop typing into them.

## Develop your regexes in regex101.com then copy them over to Flip Template to save them

Flip templates allow you to have comment lines in them using // or # character at hte start of the line. Blank lines are ignored by Flip.

The most important line to understand are the regex => substitution lines.

When you have a regex that modifies your PDF data how you like, you can paste it into the template file on a new line, and then use: " => ", that is, "space, equals, grater than, space" as the separator between the regex and whatever the substitution string for the regex will be. Be aware that additional whitespace at then end of the regex string or at the start of the substitution string (ie: directly to the left or right of the " => " separator) will affect how your regex matching works and/or what your replacement data will end up having in it, so make sure you get this right!

## What to do with your data after Flip has generated it

The data Flip generates is up to you and where the data is headed, but in the instance here of the Grun project, this data will be heading to the database(s) that the Grun Project will have. 

For the MVP of the Grun project, the next file to look at from here is called "importing_csv_data_to_MVP_UserInventory_table.md"