#!/bin/bash
CommitComment=$(zenity --entry --text "Care to comment?" --entry-text "generic new commit");
git add --all
git commit -m "$CommitComment"
git push